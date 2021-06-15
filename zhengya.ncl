begin
    filename1 = "E:/Documents/Codefield/Code_f90/zhengya/ua.dat"
    filename2 = "E:/Documents/Codefield/Code_f90/zhengya/va.dat"
    filename3 = "E:/Documents/Codefield/Code_f90/zhengya/za.dat"

    u = asciiread(filename1, (/16,20/), "float")
    v = asciiread(filename2, (/16,20/), "float")
    z = asciiread(filename3, (/16,20/), "float")
    ;printMinMax(u, 0)
    ;printMinMax(v,0)

    wks = gsn_open_wks("pdf", "E:/Documents/Codefield/Code_f90/uvz")
    gsn_define_colormap(wks, "BlueWhiteOrangeRed")

    res = True
    res@gsnDraw      = False
    res@gsnFrame     = False    
    res@gsnCenterString = "wind-ght-origin"
    resv = res

    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnLineLabelsOn = True
    res@cnLineLabelBackgroundColor = -1

    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = 5000
    res@cnMaxLevelValF = 6000
    res@cnLevelSpacingF = 50 

    resv@vcGlyphStyle = "LineArrow"
    resv@vcMinDistanceF            = 0.03
    resv@vcLineArrowHeadMaxSizeF = 0.005
    resv@vcLineArrowHeadMinSizeF = 0.005  

    resv@vcRefAnnoOn              = True
    resv@vcRefLengthF             = 0.03
    resv@vcRefAnnoBackgroundColor = "white"
    resv@vcRefAnnoPerimOn         = False
    resv@vcRefAnnoFontHeightF     = 0.015
    resv@vcRefAnnoString1On       = False
    resv@vcRefAnnoString2On       = True
    resv@vcRefAnnoSide            = "top"
    resv@vcRefAnnoOrthogonalPosF  = -0.05
    resv@vcRefAnnoParallelPosF    = 1
    resv@vcRefMagnitudeF   = 20
    resv@vcRefAnnoString2  = "20m/s"


    base = gsn_csm_contour(wks, z, res)
    plot = gsn_csm_vector(wks, u, v, resv)
    overlay(base, plot)
    draw(base)
    frame(wks) 

end