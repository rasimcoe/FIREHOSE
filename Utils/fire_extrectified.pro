pro fire_extrectified, rectflux

  flux_orders = fltarr(2048, 21)
  err_orders  = fltarr(2048, 21)
  wv_orders  = fltarr(2048, 21)

  for i=0L, 20L do begin

     if (i EQ 0) then begin
        hdr = headfits(rectflux)
     endif
     
     flux = xmrdfits(rectflux, i)
     recterr = sxpar(hdr, "ERRFIL")
     rectwv  = sxpar(hdr, "WVFIL")
     recttrc = sxpar(hdr, "TRCFIL")

     err  = xmrdfits(recterr, i)
     wv   = xmrdfits(rectwv, i)
     trc  = xmrdfits(recttrc, i)

     mask = 0.*flux + 1.0
     sz = size(flux)
     ny = sz[2]
     trace = fltarr(ny) + 15.0
     
     mask = 0.*flux + 1.0
     delta = median(wv[*,ny/2.]-wv[*,ny/2.-1])
     for j=0, ny-1 do begin
        bad = where(abs(wv[*,j]-median(wv[*,j])) GT delta, nbad)
        if(nbad GT 0) then begin
           mask[bad,j] = 0
        endif
     endfor

;     spec = long_extract_optimal(wv, flux, 1./err^2, trc, mask,
;     err^2, trace, box_rad=4.0)

     ; Fit a profile model (gaussian, coubld improve if we like) to
     ; the tennuric standard
     smash = total(trc,2)
     fit = gaussfit(findgen(n_elements(smash)),smash,coeffs, nterms=3)
     fit /= total(fit)
     profimg = (fltarr(ny) + 1.0) ## fit
     profimg[where(profimg LT 0.1 * max(fit))] = 0

     spec = fltarr(ny)
     err1d  = fltarr(ny)
     wave  = fltarr(ny)

     ; For each pixel, total up all pixels with a non-zero profile value.
     ; Note that up above, we set all pixels with a value of <10%
     ; the profile peak value to zero.  Corresponds to a +/-2 sigma 
     ; aperture, approximately.  Uses Horne extraction
     for j=0, ny-1 do begin
        spec[j] = total(profimg[*,j] * flux[*,j] / err[*,j]^2) / total(profimg[*,j]^2/err[*,j]^2)
;        wave[j] = total(profimg[*,j] * wv[*,j] / err[*,j]^2) / total(profimg[*,j]^2/err[*,j]^2)
        wave[j] = total(profimg[*,j] * wv[*,j]) / total(profimg[*,j])
        err1d[j]  = sqrt(1./total(profimg[*,j]^2/double(err[*,j])^2))
     endfor

     flux_orders[*,i] = spec
     err_orders[*,i] = err1d
     wv_orders[*,i] = wave

  endfor

  ; now rebin the spectrum onto a common wavelength grid, with inverse
  ; variance weighting to opimize SNR.  Splices the individual orders
  ; into a single continuous spectrum. 

  velpix=12.5d
  cdelt = velpix/299792.458d/alog(10.0d)
  npix=alog10(26000./8000.)/cdelt
  bigwv = 10^(alog10(8000.d) + dindgen(npix)*cdelt)
  bigfx = double(0.0 * bigwv)
  bigwt = bigfx
  biger = bigfx

  ; Perform a linear interpolation onto the big grid and then 
  ; perform an inverse variance weighted coadd.

  for i=0L, 20L do begin
     ind = sort(wv_orders[*,i])
     interpspec, wv_orders[ind,i], flux_orders[ind,i], bigwv, tmp, tmperr, yaerror=err_orders[ind,i]
     gd = where(tmperr GT 0 AND finite(tmperr) EQ 1, ngd)
     thiswt = 1./double(tmperr[gd])^2
     bigwt[gd] += thiswt
     bigfx[gd] += thiswt * tmp[gd]
     biger[gd] += 1./double(tmperr[gd])^2

     bad = where(finite(bigwt(gd)) NE 1, nbad)
     if (nbad GT 0) then stop
  endfor

  bigfx = bigfx / (bigwt + (bigwt EQ 0.0))
  biger = sqrt( biger / (bigwt^2))

  bigfx[where((bigwv GT 13150 AND bigwv LT 14100) OR $
              (bigwv GT 17600 AND bigwv LT 19300) OR (bigwv LT 8450) OR (bigwv GT 23000))] = 0 

  biger[where((bigwv GT 13150 AND bigwv LT 14100) OR $
              (bigwv GT 17600 AND bigwv LT 19300) OR (bigwv LT 8450) OR (bigwv GT 23000))] = 0 
              

  newhdr = hdr
  sxdelpar, newhdr, "NAXIS"
  sxdelpar, newhdr, "NAXIS1"
  sxdelpar, newhdr, "NAXIS2"
  sxdelpar, newhdr, 'XTENSION'
  sxdelpar, newhdr, 'PCOUNT'
  sxdelpar, newhdr, 'GCOUNT'
  sxdelpar, newhdr, 'TFIELDS'
  sxdelpar, newhdr, 'EXTEND'
  sxdelpar, newhdr, 'BITPIX'

  sxaddpar, newhdr, "CRVAL1", alog10(8000.d)
  sxaddpar, newhdr, "CDELT1", cdelt
  sxaddpar, newhdr, "CRPIX1", 1
  sxaddpar, newhdr, "CTYPE1", 'LINEAR'
  sxaddpar, newhdr, 'DC-FLAG', 1

  arr = strsplit(rectflux, '.', /extract)
  newname = 'FSpec/'+strtrim(arr[0],2)+'_F2d.fits'
  newerr = 'FSpec/'+strtrim(arr[0],2)+'_E2d.fits'
  
  mwrfits, bigfx, newname, newhdr, /create
  mwrfits, biger, newerr, newhdr, /create

  print, "Writing output to:"
  print, newname
  print, newerr
  print, ' '
  print, 'All done!'

end
