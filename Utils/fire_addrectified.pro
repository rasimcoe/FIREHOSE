;
;
; This is a supporting function that co-adds the 2D 
; rectified order spectra with a little intelligent masking
; and error handling.  Uses inverse variance weighting.
; 
function ras_2dcomb, fxarr, erarr, OUTERR=outerr

  dims = size(fxarr)
  nstack = dims[3]

  medimg = median(fxarr,dimension=3)
  sigimg = 0.0 * medimg

  ; There are sometimes bad NaNs that we need to eliminate
  bad = where(finite(erarr) EQ 0, nbad)

  ; These are the weights
  wt = 1 / erarr^2
  wtd = wt
  mask = 0.0 * erarr + 1.0

  if (nbad GT 0) then begin
     mask[bad] = 0.0
  endif

  ; Some sigma rejection, use median statistics for robustness
;
; Procedure: 1) for each row (i.e. single wavelength) calculate the
;     quartiles of the pixel values in that rown, and convert to 
;     a Gaussian sigma.  Reject pixels that are 7\sigma above the 
;     median as cosmic rays.  WARNING this could reject object pixels
;     in certain circumstances where there are high SNR.  
; 
;     2) Then, go through the stack in the image direction and throw
;        out pixels that are +/- 3 sigma from the median.
; 
  ny = dims[2]
  for i=0, ny-1 do begin
     tmp    = reform(fxarr[*,i,*])
     sorted = tmp[sort(tmp)]
     med    = sorted[n_elements(sorted) * 0.50]
     iqr    = sorted[n_elements(sorted) * 0.75]-sorted[n_elements(sorted) *0.25]
     sigma  = iqr / 1.349 ; Magic number for a gaussian dist.
     sigimg[*,i] = sigma

     for j=0, nstack-1 do begin
        hot = where(abs(fxarr[*,i,j]-med) GT 7.0 * sigma, nhot)
        if (nhot GT 0) then begin
;           fxarr[hot,i,j] = 0.0
           mask[hot,i,j] = 0.0
        endif
     endfor
  endfor

  for i=0, nstack-1 do begin
     bad = where(abs(fxarr[*,*,i]-medimg) GT 3.0 * sigimg, nbad)
     if (nbad GT 0) then begin
        mask[bad,i] = 0.0
     endif
  endfor
  
  ; Weighted mean (inverse variance)
  outflux = total(fxarr * wt * mask, 3)
  denom = total(wtd * mask, 3)
  outflux /= denom
  
  outerr = sqrt(total(erarr^2 * wt^2 * mask, 3) / (denom + (denom EQ 0.0))^2)

  bb = where(denom EQ 0 OR finite(outflux) NE 1)
  outflux[bb] = 0.0
  outerr[bb] = 0.0

  return, outflux

end

;;;
;
;
;
; TELLFILE keyword: if passed, this fluxes the output counts.  The
; filename passed to tellfile should be the the 
;
; Tellspec_0001_0002_tellspec.fits
;
; file for a corresponding telluric exposure in the Object directory.
;
pro fire_addrectified, obj_id, tellfile=tellfile, clobber=clobber

  firestrct = xmrdfits("firestrct.fits",1)

  scistrct = firestrct[where(firestrct.obj_id EQ obj_id)]

  A = scistrct[where(scistrct.posslit EQ 'A')]
  B = scistrct[where(scistrct.posslit EQ 'B')]

  Rect_a = strarr(n_elements(a))
  rect_b = strarr(n_elements(b))
  recterr_a = strarr(n_elements(a))
  recterr_b = strarr(n_elements(b))

  ; This grabs the corresponding telluric cal spectrum for flux cal
  ; Must run telluric process in Firehose to generate this.
  if (keyword_set(TELLFILE)) then begin
     tt = xmrdfits(tellfile)
     tellwv = reform(tt[*,0,*])
     tellfx = reform(tt[*,1,*])
     teller = reform(tt[*,2,*])
  endif

;;;;;;;;;;;;;;;  Rectify the A slit position exposures ;;;;;;;;;;;;

  exptime_A = 0

  for i=0, n_elements(A)-1 do begin

     arr = strsplit(A[i].fitsfile,"_.",/extract)
     number = arr[1]

     finalfile = 'Final/f'+strtrim(number,2)+'.fits.gz'
     rectfile  = 'Rectify/Rect_'+strtrim(number,2)+'.fits'
     recterr   = 'Rectify/RectErr_'+strtrim(number,2)+'.fits'
     rect_a[i] = rectfile
     recterr_a[i] = recterr
     orderfile = A[i].orderfile
     rectwv    = 'Rectify/Rectwv_'+strtrim(number,2)+'.fits'
     arcfile   = 'Arcs/ArcImg'+strtrim(A[i].arcs,2)+'.fits.gz'
     stackfile_a = strtrim(a[0].object,2)+'A.fits'
     stackwv   = strtrim(a[0].object,2)+'wv.fits'
     stacker_a = strtrim(a[0].object,2)+'Ae.fits'
     tfiles    = fire_get(A[i],/TFILES)
     slitpos_tell = strarr(n_elements(tfiles))

     ; Find the "A" position files for the telluric
     for j=0, n_elements(tfiles)-1 do begin
;        slitpos_tell[j] = sxpar(headfits(tfiles[j]), "SLITPOS")
        here = where(tfiles[j] EQ strtrim(firestrct.rawpath+firestrct.fitsfile,2))
        slitpos_tell[j] = firestrct[here].posslit
     endfor

     ; Find the tellurics with an "A" slit posn
     match = where(A[i].posslit EQ slitpos_tell, nmatch)
     if (nmatch GE 1) then begin
        arr = strsplit(tfiles[match[0]],"_.",/extract)
        number = arr[1]
        finaltell = 'Final/f'+strtrim(number,2)+'.fits.gz'
        recttell  = 'Rectify/Rect_'+strtrim(number,2)+'.fits'
     endif else begin
        print, "Warning: no telluric found for tracing this file!"
        print, " "
        print, "Typically, this is because the A/B slit posn in the raw data header is incorrect."
        print, "These are the telluric trace files for this object expsure:"
        print, tfiles
        print, "Now these are the slit positions firehose thinks these have:"
        print, slitpos_tell
        print, "If this is incorrect, you probably need to edit the firestrct_script.txt file"
        print, "and reconstitute the fire structure to fix these tags and try again."
        stop
     endelse

     print, finalfile
     print, rectfile
     print, orderfile
     print, rectwv
     print, arcfile
     print, stackfile_a
     print, finaltell

     exptime_A += sxpar(headfits(finalfile),"EXPTIME")

     if (file_test(rectfile) EQ 1 and (not keyword_set(CLOBBER))) then begin
        continue
     endif else begin
        ; Rectify the sky-subtracted science frame
        fire_rectify, finalfile, orderfile, rectfile, $
                      /arcsec, /skysub, waveoutpt=rectwv, $
                      waveimg=arcfile, outerr=recterr
           ; Rectify the non-sky-subtracted telluric (for tracing purposes)
        fire_rectify, finaltell, orderfile, recttell, /arcsec
     endelse

  endfor

  recttell_a = recttell

;;;;;;;;;;;;;;;  Rectify the B slit position exposures ;;;;;;;;;;;;

  exptime_B = 0

  for i=0, n_elements(B)-1 do begin

     arr = strsplit(B[i].fitsfile,"_.",/extract)
     number = arr[1]

     finalfile = 'Final/f'+strtrim(number,2)+'.fits.gz'
     rectfile  = 'Rectify/Rect_'+strtrim(number,2)+'.fits'
     recterr   = 'Rectify/RectErr_'+strtrim(number,2)+'.fits'
     rect_b[i] = rectfile
     recterr_b[i] = recterr
     orderfile = B[i].orderfile
     rectwv    = 'Rectify/Rectwv_'+strtrim(number,2)+'.fits'
     arcfile   = 'Arcs/ArcImg'+strtrim(B[i].arcs,2)+'.fits.gz'
     stackfile_b = strtrim(b[0].object,2)+'B.fits'
     stacker_b   = strtrim(a[0].object,2)+'Be.fits'
;     tfiles    = fire_get(B[i],/TFILES)
     tfiles    = strtrim(B[i].rawpath,2)+strtrim(B[i].objtracefile,2)
     slitpos_tell = strarr(n_elements(tfiles))
     for j=0, n_elements(tfiles)-1 do begin
;        slitpos_tell[j] = sxpar(headfits(tfiles[j]), "SLITPOS")
        here = where(tfiles[j] EQ strtrim(firestrct.rawpath+firestrct.fitsfile,2))
        slitpos_tell[j] = firestrct[here].posslit
     endfor
     match = where(B[i].posslit EQ slitpos_tell, nmatch)
     if (nmatch GE 1) then begin
        arr = strsplit(tfiles[match[0]],"_.",/extract)
        number = arr[1]
        finaltell = 'Final/f'+strtrim(number,2)+'.fits.gz'
        recttell  = 'Rectify/Rect_'+strtrim(number,2)+'.fits'
     endif else begin
        print, "Warning: no telluric found for tracing this file!"
        stop
     endelse


     print, finalfile
     print, rectfile
     print, orderfile
     print, rectwv
     print, arcfile
     print, stackfile_b

     exptime_B += sxpar(headfits(finalfile),"EXPTIME")

     if (file_test(rectfile) EQ 1 and (not keyword_set(CLOBBER))) then begin
        continue
     endif else begin
        fire_rectify, finalfile, orderfile, rectfile, $
                      /arcsec, /skysub, waveoutpt=rectwv, $
                      waveimg=arcfile, outerr=recterr
        ; Rectify the non-sky-subtracted telluric (for tracing purposes)
;        fire_rectify, finaltell, orderfile, recttell, /arcsec
     endelse

  endfor

  recttell_b = recttell

 ;;;;;;;; Now, combine the A and B rectified exposures ;;;;;;;;;;;

  ; Must go order by order
  ; Stack up and run djs_avsigclip

  for i=0L, 20 do begin

     for j=0, n_elements(A)-1 do begin
        rectfile = rect_a[j]
        recterr = recterr_a[j]
        if (j EQ 0) then begin
           hdr_a = headfits(rectfile)
           buf = xmrdfits(rectfile,i)
           exptime = sxpar(hdr_a, "EXPTIME")
           buf /= exptime
           sz = size(buf)
           arr = fltarr(sz[1],sz[2],n_elements(A))
           err_arr = fltarr(sz[1],sz[2],n_elements(A))
           arr[*,*,j] = buf
           err_arr[*,*,j] = xmrdfits(recterr,i) / exptime
        endif else begin
           buf = xmrdfits(rectfile,i)
           hh = headfits(rectfile)
           exptime = sxpar(hh, "EXPTIME")
           buf /= exptime
           arr[*,*,j] = buf
           err_arr[*,*,j] = xmrdfits(recterr,i) / exptime
        endelse
     endfor
     
     if (n_elements(A) GT 1) then begin
        Acomp = ras_2dcomb(arr, err_arr, outerr=Aerr)
     endif else begin
        Acomp = arr[*,*,0]
        Aerr = err_arr[*,*,0]
     endelse
     

     for j=0, n_elements(B)-1 do begin
        rectfile = rect_b[j]
        if (j EQ 0) then begin
           hdr_b = headfits(rectfile)
           buf = xmrdfits(rectfile,i)
           exptime = sxpar(hdr_b, "EXPTIME")
           buf /= exptime
           sz = size(buf)
           arr = fltarr(sz[1],sz[2],n_elements(B))
           err_arr = fltarr(sz[1],sz[2],n_elements(B))
           arr[*,*,j] = buf
           err_arr[*,*,j] = xmrdfits(recterr,i) / exptime
        endif else begin
           buf = xmrdfits(rectfile,i)
           hh = headfits(rectfile)
           exptime = sxpar(hh, "EXPTIME")
           buf /= exptime
           arr[*,*,j] = buf
           err_arr[*,*,j] = xmrdfits(recterr,i) / exptime
        endelse
     endfor

     if (n_elements(B) GT 1) then begin
        Bcomp = ras_2dcomb(arr, err_arr, outerr=Berr)
     endif else begin
        Bcomp = arr[*,*,0]
        Berr = err_arr[*,*,0]
     endelse

     wv = xmrdfits(rectwv,i)

     tmp = i

     ;; Now flux-calibrate the output if so desired (TELLFILE passed in)
     ;; Aout and Bout are the co-added frames in cts/sec
     Aout = AComp
     Bout = BComp
     if (keyword_set(TELLFILE)) then begin

        ; The tellspec outputs of xtellcor are in units of
        ; erg/cm2/sec/A/DN (see headers in output of tellspec files)
        ww = tellwv[*,i]
        ff = tellfx[*,i]

        w1d = reform(wv[sz[1]/3.,*])
        interpspec, ww*10000., ff, w1d, f1d
        sensarr = f1d ## (fltarr(sz[1])+1.)
        Aout *= sensarr
        Bout *= sensarr
        Aerr *= sensarr
        Berr *= sensarr
     endif

     sxaddpar, hdr_a, "EXPTIME", exptime_A
     sxaddpar, hdr_b, "EXPTIME", exptime_B
     sxaddpar, hdr_a, "TRCFIL", recttell_a
     sxaddpar, hdr_b, "TRCFIL", recttell_b
     sxaddpar, hdr_a, "ERRFIL", stacker_a
     sxaddpar, hdr_b, "ERRFIL", stacker_b
     sxaddpar, hdr_a, "WVFIL", stackwv
     sxaddpar, hdr_b, "WVFIL", stackwv

     ; Write out the results
     if (tmp EQ 0) then begin
        mwrfits, Aout, stackfile_a,hdr_a, /create
        mwrfits, Aerr, stacker_a,hdr_a, /create
        mwrfits, Bout, stackfile_b,hdr_b, /create
        mwrfits, Berr, stacker_b,hdr_b, /create
        mwrfits, wv, stackwv, /create
     endif else begin
        mwrfits, Aout, stackfile_a
        mwrfits, Aerr, stacker_a
        mwrfits, Bout, stackfile_b
        mwrfits, Berr, stacker_b
        mwrfits, wv, stackwv
     endelse

  endfor

  ; All done!

end
