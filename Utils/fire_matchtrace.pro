; Given an input fire structure, find all SCIENCE exposures and then
; determine the best OBJTRACEFILE to be used for deep exposures where
; there is no continuum expected.  Here, we make use of the TFILES
; kepword to find the nearest-time telluric matched with that object
; frame, and then cross that with the POSSLIT keyword to make sure
; that the trace for objects at the A position use telluric at A and
; vice versa.

PRO fire_matchtrace, fire, SAVE=save

  for i=0, n_elements(fire)-1 do begin

     if (fire[i].exptype EQ 'SCIENCE') then begin
        tfiles = fire_get(fire[i], /TFILES)
        slitpos_tell = strarr(n_elements(tfiles))
        for j=0, n_elements(tfiles)-1 do begin
           slitpos_tell[j] = sxpar(headfits(tfiles[j]), "SLITPOS")
        endfor

        match = where(fire[i].posslit EQ slitpos_tell, nmatch)
        if (nmatch GT 0) then begin
           arr = strsplit(tfiles[match[0]],'/', /extract)
           objtrace = arr[n_elements(arr)-1]
           fire[i].objtracefile = objtrace
           print, fire[i].object, fire[i].fitsfile,  fire[i].posslit, objtrace
        endif

     endif

  endfor

  if (keyword_set(SAVE)) then begin
     mwrfits, fire, 'firestrct.fits', /create
  endif

end
