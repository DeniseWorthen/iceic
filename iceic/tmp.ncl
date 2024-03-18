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
  ;minlat = 89.8
  ;minlat = 88.0

   ;vmag = 0.80
   ;vthn = 0.01
   vthn = 0.008

   ;rad2deg=1.0
   rad2deg=180.d0/acos(-1.0d0)

;************************************************
; read in netCDF file
;************************************************

  RT2 = "/scratch1/NCEPDEV/stmp2/Denise.Worthen/"

  cdate = "2011-10-03-00000"
  rtsrc = "tripole/"
  plotlabs = (/"develop", "fix"/)
       rts = (/"base.025.dev/", "base.025.fixno360/"/)
   icehist = "iceh_06h."+cdate+".nc"
  icefiles = addfiles(RT2+rtsrc+rts+"history/"+icehist,"r")
  print(icefiles)
  ListSetType(icefiles,"join")

  icelonname = "ULON"
  icelatname = "ULAT"
   anglename = "ANGLE"
   plottitle = "cice5"

      ut = icefiles[:]->uvel_h(:,0,:,:)
      vt = icefiles[:]->vvel_h(:,0,:,:)
  ut@lon2d=icefiles[:]->$icelonname$
  ut@lat2d=icefiles[:]->$icelatname$
  vt@lon2d=icefiles[:]->$icelonname$
  vt@lat2d=icefiles[:]->$icelatname$

  ;avar = icefiles[:]->aice_h(:,0,:,:)

  ; copy metadata, rotang here is same sign as internal CICE (-1*sin,+1*cos mom)
  rotang = icefiles[:]->$anglename$
  ;printVarSummary(rotang)
  ;printVarSummary(ut)
  ;urot = where(avar .ge. 0.8, ut, default_fillvalue(typeof(ut)))
  ;vrot = where(avar .ge. 0.8, vt, default_fillvalue(typeof(ut)))
  urot = ut
  vrot = vt
  urot = cos(rotang)*ut - sin(rotang)*vt
  vrot = sin(rotang)*ut + cos(rotang)*vt

  ;printVarSummary(urot)
;if(1.eq.0)then
  urot = where(urot@lat2d .ge. minlat , urot, default_fillvalue(typeof(ut)))
  vrot = where(vrot@lat2d .ge. minlat , vrot, default_fillvalue(typeof(vt)))
  urot@_FillValue = default_fillvalue(typeof(ut))
  vrot@_FillValue = default_fillvalue(typeof(vt))
  print("urot min = "+min(urot)+" urot max = "+max(urot))
  print("vrot min = "+min(vrot)+" vrot max = "+max(vrot))


  ut = urot
  vt = vrot
  printVarSummary(ut)

  ;udiff = urot(1,:,:) - urot(0,:,:)
  ;vdiff = vrot(1,:,:) - vrot(0,:,:)
  ;printVarSummary(udiff)

  ;udiff@lon2d=icefiles[1]->$icelonname$
  ;udiff@lat2d=icefiles[1]->$icelatname$
  ;vdiff@lon2d=icefiles[1]->$icelonname$
  ;vdiff@lat2d=icefiles[1]->$icelatname$

  ; copy metadata
  ;udiff = urot(0,:,:)
  ;vdiff = vrot(0,:,:)
  ;udiff = urot(1,:,:) - urot(0,:,:)
  ;vdiff = vrot(1,:,:) - vrot(0,:,:)
  ;ut@lat2d = icefiles[1]->$icelonname$
  ;vt@lat2d = icefiles[1]->$icelonname$
  ;ut@lon2d = icefiles[1]->$icelatname$
  ;vt@lon2d = icefiles[1]->$icelatname$
  ;printVarSummary(udiff)

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

  plot = new(1,graphic)
  oplot = new(1,graphic)

  res                      = True        ; plot mods desired
  res@gsnDraw              = False       ; do not draw until the end
  res@gsnFrame             = False       ; do not automatically advance 'frame'

  res@cnFillMode = "RasterFill"
  res@cnFillOn          = True                  ; color fill
  res@cnLinesOn         = False                 ; no contour lines
  res@cnLineLabelsOn    = False
  res@lbLabelBarOn      = False

  ;res@vcRefMagnitudeF       = vmag              ; vec magnitude
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

  uvar= ut(0,:,:)
  vvar= vt(0,:,:)
  ;uvar@lon2d = icefiles[1]->$icelonname$
  ;vvar@lon2d = icefiles[1]->$icelonname$
  ;uvar@lat2d = icefiles[1]->$icelatname$
  ;vvar@lat2d = icefiles[1]->$icelatname$

  printVarSummary(uvar)
  printVarSummary(vvar)

  ;plot = gsn_csm_vector_map_polar(wks,uvar,vvar,res)    ; create the plot

  ;plot = gsn_csm_contour_map_polar(wks,vvar,res)    ; create the plot
  ;dum(0)= gsn_add_polyline(wks,plot(0),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
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
;end if
exit
end