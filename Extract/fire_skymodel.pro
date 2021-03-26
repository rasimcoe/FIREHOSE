function fire_skymodel, SCIIMG=sciimg, SCIIVAR=sciivar, $
                        PIXIMAGE=piximage, HDR=hdr, TSET_SLITS=tset_slits, $
                        STD=std, DATA=data, VERBOSE=verbose, $
                        CHECK_INPUTS=check_inputs, shiftpix=shiftpix, CHK=chk
  
  prog_name = 'fire_skymodel.pro'
  print, ''
  print, prog_name + ": Commencing fire_skymodel..."
  
  sciimgfile = "Not read in: fed via keyword."
  piximgfile = "Not read in: fed via keyword."
  tset_slitsfile = "Not read in: fed via keyword."
  
  if keyword_set(DATA) then begin
     
     ;; if the sciimg (or sciivar) argument is not present, we must determine it
     if NOT keyword_set(sciimg) OR NOT keyword_set(sciivar) then begin
        print, prog_name, ": Generating sciimg, sciivar, and hdr..."
        sciimg = fire_get(data, /SCIIMG, STD=std, HDR=hdr, SCIIVAR=sciivar, VERBOSE=verbose, FILEREAD=sciimgfile)
     endif
     if NOT keyword_set(piximg) then begin
        print, prog_name, ": Generating piximg..."
        piximg = fire_get(data, /PIXIMAGE, VERBOSE=verbose, FILEREAD=piximgfile)
     endif
     if NOT keyword_set(TSET_SLITS) then begin
        print, prog_name, ": Generating tset_slits..."
        tset_slits = fire_get(data,/TSET_SLITS, VERBOSE=verbose, FILEREAD=tset_slitsfile)
     endif
     
  endif
  
  if NOT keyword_set(DATA) AND $
     ( NOT arg_present(sciimg) OR NOT keyword_set(PIXIMAGE) OR NOT keyword_set(TSET_SLITS) ) then begin
     fire_siren, prog_name + ": invalid input to fire_skymodel!  If keyword DATA is not passed, then the argument sciimg must be present, and the keywords PIXIMG and TSET_SLITS must be passed!  Exiting with error!"
     RETURN, 0
  endif
  
  if keyword_set(VERBOSE) OR keyword_set(CHECK_INPUTS) then begin
     print, ''
     print, "++++++++++++++++++++++++++++++++++++"
     print, prog_name, ": /VERBOSE flag passed.  Set values:"
     print, "Ran fire_proc on file: ", sciimgfile
     print, "Read in pixel image from file: ", piximgfile
     print, "Read in tset_slits from file: ", tset_slitsfile
     print, "++++++++++++++++++++++++++++++++++++"
     print, ''
     if keyword_set(CHECK_INPUTS) then RETURN, 0
  endif
  
  bsp   = 1.0                    ;tweaking sky model
  FWHM  = 2.5
  npoly = 0
  objstruct1 = long_objfind(sciimg, tset_slits = tset_slits, invvar = sciivar $
                            , skymask = skymask, objmask = objmask $
                            , nperslit = 1L, peakthresh = reduxthresh $
                            , fwhm = FWHM, ISLIT = ISLIT)
  slitmask = long_slits2mask(tset_slits)
  ;; aggresive edge masking because of order flexure
  ximg = long_slits2x(tset_slits, edgmask = edgmask,TOL_EDG=0)
  
  ; The following code is a hack to mask other objects in the
  ; slit out by hand for the purpose of sky subtraction.  User should 
  ; set the minimum and maximum range of pixels in the slit to exclude
  ; from the sky estimate.  The units are fraction of a slit, where 
  ; the left hand edge of the slit is zero, and the right edge is 1.
  ; So, the masking min and max should be fractions.
  tweak_skymask = 0
  if (tweak_skymask EQ 1) then begin
     min_mask = -1
     max_mask = -1
     stop
     if (min_mask GT 0.0 AND min_mask LT 1.0 AND $
         max_mask GT 0.0 AND max_mask LT 1.0 AND min_mask LT max_mask) then begin
        skymask[where(ximg GT min_mask AND ximg LT max_mask)] = 0
        print, 'Displaying pixels which will be used to estimate the sky for sky subtraction purposes.'
        xatv, sciimg * skymask, /block
     endif
  endif

  if (NOT keyword_set(SHIFTPIX)) then begin
     skyimage = long_skysub(sciimg, sciivar, piximage, slitmask, skymask, edgmask $
                            , bsp = bsp, ISLIT = ISLIT, npoly=npoly)
  endif else begin
     skyimage = fire_skyshift(sciimg, sciivar, piximage, slitmask, skymask, edgmask $
                              , bsp = bsp, ISLIT = ISLIT, npoly=npoly, shiftpix=shiftpix)
  endelse

;  skyimage = long_skysub(sciimg, sciivar, piximage, slitmask, skymask, edgmask $
;                         , bsp = 1, ISLIT = ISLIT, npoly=npoly)

  if (keyword_set(CHK)) then begin
     FOR iorder=0, 20 do begin
        order = iorder
        inord = where(slitmask EQ order)
        xx  = piximage[inord]
        sp  = sciimg[inord]
        mdl = skyimage[inord]
        colors=getcolor(/load)
        xxx = xx[sort(xx)]
        mdl2 = mdl[sort(xx)]
        x_splot, xx, sp, xtwo=xxx, ytwo=mdl2, /block, psym1=3
     ENDFOR
  endif

  print, prog_name + ": ...completed fire_skymodel."
  print, ''
  return, skyimage
  
end
