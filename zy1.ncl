begin
    filename1 = "E:/Documents/Codefield/Code_f90/zhengya/uc1.dat"
    filename2 = "E:/Documents/Codefield/Code_f90/zhengya/vc1.dat"
    filename3 = "E:/Documents/Codefield/Code_f90/zhengya/zc1.dat"
    filename4 = "E:/Documents/Codefield/Code_f90/zhengya/uc2.dat"
    filename5 = "E:/Documents/Codefield/Code_f90/zhengya/vc2.dat"
    filename6 = "E:/Documents/Codefield/Code_f90/zhengya/zc2.dat"

    u1 = asciiread(filename1, (/16,20/), "float")
    v1 = asciiread(filename2, (/16,20/), "float")
    z1 = asciiread(filename3, (/16,20/), "float")
    u2 = asciiread(filename4, (/16,20/), "float")
    v2 = asciiread(filename5, (/16,20/), "float")
    z2 = asciiread(filename6, (/16,20/), "float")

    u3 = u2-u1
    v3 = v2-v1
    z3 = z2-z1
    printMinMax(u1, 0)
    printMinMax(v1,0)
    printMinMax(z1,0)
    printMinMax(u2, 0)
    printMinMax(v2,0)
    printMinMax(z2,0)
    printMinMax(u3,0)
    printMinMax(v3,0)
    printMinMax(z3,0)

    wks = gsn_open_wks("pdf", "E:/Documents/Codefield/Code_f90/uvz1")
    gsn_define_colormap(wks, "BlueWhiteOrangeRed")

    res = True
    res@gsnDraw      = False
    res@gsnFrame     = False    
    res@gsnCenterString = "wind-ght-SJPH"
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
    ;res2 = res
    ;resv2 = resv
    ;res2@cnMinLevelValF = 2500
    ;res2@cnMaxLevelValF = 9000
    ;res2@cnLevelSpacingF = 325
    ;resv2@vcRefMagnitudeF   = 100
    ;resv2@vcRefAnnoString2  = "100m/s"
    res3 = res
    resv3 = resv
    res3@cnMinLevelValF = -500
    res3@cnMaxLevelValF = 500
    res3@cnLevelSpacingF = 50
    resv3@vcRefMagnitudeF   = 10
    resv3@vcRefAnnoString2  = "10m/s"

    base1 = gsn_csm_contour(wks, z1, res)
    plot1 = gsn_csm_vector(wks, u1, v1, resv)
    overlay(base1, plot1)
    base2 = gsn_csm_contour(wks, z2, res)
    plot2 = gsn_csm_vector(wks, u2, v2, resv)
    overlay(base2, plot2)
    base3 = gsn_csm_contour(wks, z3, res3)
    plot3 = gsn_csm_vector(wks, u3, v3, resv3)
    overlay(base3, plot3)

    draw(base1)
    frame(wks)
    draw(base2)
    frame(wks)
    draw(base3)
    frame(wks) 

end