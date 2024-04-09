;************************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;************************************************
;
;************************************************

  batch = "no"

  hemi="NH"
  ;hemi="SH"

  ;minlat = 66.0
  ;minlat = 76.0
  ;minlat = 89.8
  minlat = 88.0

   vmag = 0.50
   vthn = 0.01
   ;vthn = 0.0

   rad2deg=1.0
   ;rad2deg=180.d0/acos(-1.0d0)
begin
;************************************************
; read in netCDF file
;************************************************

  RT1 = "/scratch3/NCEPDEV/stmp1/Denise.Worthen/"
  RT2 = "/scratch3/NCEPDEV/stmp2/Denise.Worthen/"
  RT3 = "/scratch3/NCEPDEV/stmp3/Denise.Worthen/"
  RT4 = "/scratch4/NCEPDEV/stmp4/Denise.Worthen/"

  ;cdate = "2013040100"
  ;cdate = "2011090100"
  cdate = "2011-10-01-43200"

  ;icdir = "/scratch2/NAGAPE/epic/UFS-WM_RT/NEMSfv3gfs/BM_IC-20220207/2013040100/cpc/"
  ;cice5 = "cice5_model_0.25.res_"+cdate+".nc"

  icdir="/scratch1/NCEPDEV/nems/Denise.Worthen/WORK/icgit/"
  ;cice5="testdstCt_ew.nc"
  cice5="testdstBu_ij.nc"
  ;cice5="testdstBu_ew.nc"
  ;icdir = "/scratch2/NAGAPE/epic/UFS-WM_RT/NEMSfv3gfs/develop-20240305/datm_cdeps_mx025_gefs_intel/RESTART/"
  ;cice5 = "iced."+cdate+".nc"
  ; this copy contains the center lats and lons in degrees, which
  ; is more convenient for plotting
  dirsrc = "/scratch1/NCEPDEV/climate/Denise.Worthen/grids-20240311/"
  ; load the grid file
  ;gr = addfile(dirsrc+"grid_cice_NEMS_mx025.nc","r")
  ;gr = addfile(dirsrc+"grid_cice_NEMS_mx100.nc","r")
  gr = addfile(dirsrc+"tripole.mx100.nc","r")

  icelonname = "lonBu"
  icelatname = "latBu"
   anglename = "angle"
   plottitle = cice5

  ;print(icdir+cice5)
  ifile = addfile(icdir+cice5,"r")
  wkstitle="restarts"+cdate

  rotang = gr->$anglename$
      ut = ifile->vecdst(0,:,:)
      vt = ifile->vecdst(1,:,:)
  ;ut = ifile->uvel
  ;vt = ifile->vvel
      ;printVarSummary(ut)
      ;print(ut(1079,:)+"   "+vt(1079,:))

  ;ut@_FillValue =  default_fillvalue("double")
  ;vt@_FillValue =  default_fillvalue("double")

  ut@lon2d=gr->$icelonname$*rad2deg
  ut@lat2d=gr->$icelatname$*rad2deg
  vt@lon2d=gr->$icelonname$*rad2deg
  vt@lat2d=gr->$icelatname$*rad2deg

  ; copy metadata, rotang here is same sign as internal CICE (-1*sin,+1*cos mom)
  urot = ut
  vrot = vt
  ;urot = 0.0
  ;vrot = 0.0
  urot =  cos(rotang)*ut - sin(rotang)*vt
  vrot =  sin(rotang)*ut + cos(rotang)*vt

  ;printVarSummary(urot)

  urot = where(urot@lat2d .ge. minlat, urot, 0.0)
  vrot = where(vrot@lat2d .ge. minlat, vrot, 0.0)
  print("urot min = "+min(urot)+" urot max = "+max(urot))
  print("vrot min = "+min(vrot)+" vrot max = "+max(vrot))

  levs = (/-0.2, 0.2, .02/)

;************************************************
; add a poly line showing the polar seem along
; great circle path from (lon1,lat1) to (lon2,lat2)
;************************************************

   lon1=60.00
   lon2=-120.00
   lat1=60
   lat2=60

;************************************************
; define the map
;************************************************
  wkstype="X11"
  ;wkstype="png"
  wkstype@wkWidth=800
  wkstype@wkHeight=800
  wks   = gsn_open_wks (wkstype,"cice_uv")              ; send graphics to X11 window

  plot = new(2,graphic)
  ;dum = new(2,graphic)

  res                      = True        ; plot mods desired
  res@gsnDraw              = False       ; do not draw until the end
  res@gsnFrame             = False       ; do not automatically advance 'frame'

  res@cnFillOn          = True                  ; color fill
  res@cnLinesOn         = False                 ; no contour lines
  res@cnLineLabelsOn    = False
  res@lbLabelBarOn      = True

  res@vcRefMagnitudeF       = vmag              ; vec magnitude
  res@vcRefLengthF          = 0.10              ; size of ref vector
  ;res@vcMonoLineArrowColor  = True             ; multiple colors desired
  res@vcMinDistanceF        = vthn              ; thins arrows near pole
  ;res@vcLineArrowThicknessF = 1.5
  res@vcGlyphStyle          = "CurlyVector"
  ;res@vcPositionMode = "ArrowTail"

  res@gsnPolar   =  hemi                         ; specify the hemisphere
  if(hemi .eq. "NH")then
   res@mpMinLatF    =  minlat                          ; minimum lat to plot
   res@mpCenterLonF = -45
  end if
  if(hemi .eq. "SH")then
   res@mpMaxLatF    = -50                            ; minimum lat to plot
   res@mpCenterLonF = 0.0
  end if

  res@cnLevelSelectionMode = "ManualLevels"     ; set the contour levels with the following 3 resources
  res@cnMinLevelValF  = levs(0)                      ; set the minimum contour level
  res@cnMaxLevelValF  = levs(1)                      ; set the maximum contour level
  res@cnLevelSpacingF = levs(2)                      ; set the interval between contours
  res@cnFillPalette        = "BlueWhiteOrangeRed"              ; choose colormap

;************************************************
; Resources for the polyline
;************************************************

  res@mpGreatCircleLinesOn = True
  pres                  = True           ; polyline mods desired
  pres@gsLineThicknessF = 2.0            ; line thickness
  pres@gsLineColor      = "black"         ; color of lines
  pres@gsLineDashPattern =  1            ; dashes
  ;pres@tfPolyDrawOrder = "postDraw"

; ***************************************************
; panel resources
; **************************************************

  resP                       = True                 ; modify the panel plot
  resP@gsnMaximize           = True                 ; maximize in frame
  resP@gsnFrame              = False                ; don't advance frame yet
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelLabelBar      = False                 ; turn on a common labelbar for the entire panel plot

  ;resP@gsnPanelMainFontHeightF = 0.015           ; new resource added in NCL V6.4.0
  resP@gsnPanelMainString = plottitle

;************************************************
;
;************************************************

  ;plot = gsn_csm_vector_scalar_map_polar(wks,urot,vrot,avar,res)    ; create the plot
  ;plot = gsn_csm_vector_map_polar(wks,urot,vrot,res)    ; create the plot
  plot(0) = gsn_csm_contour_map_polar(wks,urot,res)    ; create the plot
  plot(1) = gsn_csm_contour_map_polar(wks,vrot,res)
  ;dum(0) gsn_add_polyline(wks,plot(0),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
  ;dum(1) gsn_add_polyline(wks,plot(1),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
  ;draw(plot)
  ;frame(wks)

  ;plot = gsn_csm_contour_map_polar(wks,vvar,res)    ; create the plot
  dum= gsn_add_polyline(wks,plot(0),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
  dum= gsn_add_polyline(wks,plot(1),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
  ;draw(plot)
  ;frame(wks)

  ;plot = gsn_csm_contour_map_polar(wks,avar,res)    ; create the plot
  ;dum= gsn_add_polyline(wks,plot,(/lon1, lon2/) , (/lat1, lat2/) ,pres)

;************************************************
; make the panel
;************************************************

  gsn_panel(wks,plot,(/1,2/),resP)

  frame(wks)

exit
end