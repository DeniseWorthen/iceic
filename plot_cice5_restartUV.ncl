;************************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;************************************************
;
;************************************************
begin

  hemi="NH"
  ;hemi="SH"

  minlat = 66.0
  ;minlat = 76.0
  ;minlat = 88.0
  ;minlat=80.0
  ;vmag = 0.5
  ;vthn = 0.015
  vthn = 0.008
  vmag = 1.0

  rad2deg=1.0
  ;rad2deg=180.d0/acos(-1.0d0)

  plot_src = 0

;************************************************
; read in netCDF file
;************************************************

  RT1 = "/scratch3/NCEPDEV/stmp1/Denise.Worthen/"
  RT2 = "/scratch3/NCEPDEV/stmp2/Denise.Worthen/"
  RT3 = "/scratch3/NCEPDEV/stmp3/Denise.Worthen/"
  RT4 = "/scratch4/NCEPDEV/stmp4/Denise.Worthen/"
  ; ics from baseline
  ICBL = "/scratch2/NAGAPE/epic/UFS-WM_RT/NEMSfv3gfs/input-data-20221101/CICE_IC/"

  cres = "025"
  icelonname = "lonBu"
  icelatname = "latBu"
  anglename = "angle"
  plottitle = "restart"

  if (plot_src .eq. 1)then
    grdsrc="/scratch1/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data/"+cres+"/"
    gr = addfile(grdsrc+"tripole.mx"+cres+".nc","r")

    ;icdir="/scratch1/NCEPDEV/stmp2/Denise.Worthen/tripole/gefs.rundir/"
    ;fname="cice_model.res.nc"
    ;ifile = addfile(icdir+fname,"r")
    ;print(ifile)
    ;vmag = 1.0
    icdir="/scratch1/NCEPDEV/stmp2/Denise.Worthen/tripole/fix.mx025/"
    fname = "iced.2011-10-03-00000.nc"
    ifile = addfile(icdir+"RESTART/"+fname,"r")
    vmag = 1.0

      wkstitle=fname
      ut = ifile->uvel(:,:)
      vt = ifile->vvel(:,:)
    end if

    if (plot_src .eq. 0)then
      icdir="./"
      ifile = addfile(icdir+"testdstCt_ij.nc","r")
      grdsrc="/scratch1/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data/100/"
      gr = addfile(grdsrc+"tripole.mx100.nc","r")

      ut = ifile->vecdst(0,:,:)
      vt = ifile->vecdst(1,:,:)
    end if
    print("ut min = "+min(ut)+" ut max = "+max(ut))
    print("vt min = "+min(vt)+" vt max = "+max(vt))

    ut@lon2d=gr->$icelonname$
    ut@lat2d=gr->$icelatname$
    vt@lon2d=gr->$icelonname$
    vt@lat2d=gr->$icelatname$

    ;maskvar = gr->wet(:,:)
    maskvar = todouble(gr->wet(:,:))
    maskvar@lon2d = gr->lonBu
    maskvar@lat2d = gr->latBu
    printVarSummary(maskvar)
 ; copy metadata, rotang here is same sign as internal CICE (-1*sin,+1*cos mom)
    rotang = gr->$anglename$
    printVarSummary(rotang)
    printVarSummary(ut)
    urot = ut
    vrot = vt
    urot = cos(rotang)*ut - sin(rotang)*vt
    vrot = sin(rotang)*ut + cos(rotang)*vt

    urot = where(urot@lat2d .ge. minlat .and. maskvar .eq. 1.0, urot, 0.0)
    vrot = where(vrot@lat2d .ge. minlat .and. maskvar .eq. 1.0, vrot, 0.0)
    print("urot min = "+min(urot)+" urot max = "+max(urot))
    print("vrot min = "+min(vrot)+" vrot max = "+max(vrot))

   ;urot = where(maskvar .eq. 0.0, default_fillvalue(typeof(urot)), urot)
   ;vrot = where(maskvar .eq. 0.0, default_fillvalue(typeof(vrot)), vrot)


    levs = (/-0.2, 0.2, .02/)
  ;levs = (/-0.05, 0.05, .001/)


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
  dum = new(2,graphic)

  res                      = True        ; plot mods desired
  res@gsnDraw              = False       ; do not draw until the end
  res@gsnFrame             = False       ; do not automatically advance 'frame'

  res@cnFillMode = "RasterFill"
  res@cnFillOn          = True                  ; color fill
  res@cnLinesOn         = False                 ; no contour lines
  res@cnLineLabelsOn    = False
  res@lbLabelBarOn      = False

  res@vcRefMagnitudeF       = vmag              ; vec magnitude
  res@vcRefLengthF          = 0.10              ; size of ref vector
  ;res@vcMonoLineArrowColor  = True             ; multiple colors desired
  res@vcMinDistanceF        = vthn              ; thins arrows near pole
  ;res@vcLineArrowThicknessF = 1.5
  ;res@vcGlyphStyle          = "CurlyVector"
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

  ;res@cnMissingValPerimOn     = False            ; turn on the missing value perimeter
  res@cnMissingValFillColor  = "gray"

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
  resP@gsnPanelLabelBar      = True                 ; turn on a common labelbar for the entire panel plot

  ;resP@gsnPanelMainFontHeightF = 0.015           ; new resource added in NCL V6.4.0
  resP@gsnPanelMainString = plottitle

;************************************************
;
;************************************************

  ;plot = gsn_csm_vector_scalar_map_polar(wks,urot,vrot,avar,res)    ; create the plot
  plot(0) = gsn_csm_vector_map_polar(wks,urot,vrot,res)    ; create the plot
  ;plot(0) = gsn_csm_contour_map_polar(wks,urot,res)    ; create the plot
  ;plot(1) = gsn_csm_contour_map_polar(wks,vrot,res)
  ;dum(0) gsn_add_polyline(wks,plot(0),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
  ;dum(1) gsn_add_polyline(wks,plot(1),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
  ;draw(plot)
  ;frame(wks)

  ;plot = gsn_csm_contour_map_polar(wks,vvar,res)    ; create the plot
  dum(0)= gsn_add_polyline(wks,plot(0),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
  ;dum(1)= gsn_add_polyline(wks,plot(1),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
  ;draw(plot)
  ;frame(wks)

  ;plot = gsn_csm_contour_map_polar(wks,avar,res)    ; create the plot
  ;dum= gsn_add_polyline(wks,plot,(/lon1, lon2/) , (/lat1, lat2/) ,pres)

;************************************************
; make the panel
;************************************************

  ;gsn_panel(wks,plot,(/1,2/),resP)
  gsn_panel(wks,plot,(/1,1/),resP)
  frame(wks)

exit
end