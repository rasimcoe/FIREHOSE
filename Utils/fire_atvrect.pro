pro fire_atvrect, filenm, extension

  A = transpose(xmrdfits(filenm, extension))
  B = transpose(xmrdfits(filenm, extension))
  wv = transpose(xmrdfits(filenm, extension))

  xatv, A-B, wv=wv, /block

;  profimg = A * 0
;  sz = size(profimg)
;  fwhm = 5.0
;  sigma = fwhm / 2.35
;  i0 = sz[2]*0.25
;  for i=0, sz[2]-1 do begin
;     profimg[*,i] = 1./(sqrt(2.*!PI)*sigma) * exp(-1.0 * (i-i0)^2/(2.*sigma^2))
;  endfor

;  spec = total(A*profimg,2)
;  x_splot, spec, /block

end
