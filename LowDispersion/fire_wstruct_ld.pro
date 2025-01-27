;+
; NAME:
;   long_wstruct
;
; PURPOSE:
;   Create a structure which guides the data reduction process.  Parts
;   of the structure are initialized according to the instrument and
;   its configuration.
;
; CALLING SEQUENCE:
;
; INPUTS:
;   hdr  -- Image header
;
; OPTIONAL INPUTS:
;                
; OUTPUTS:
;  Returns a structure describing the instrument configuration which
;  is used to guide the reduction steps.
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Mar-2005  Written by S. Burles (MIT), David Schlegel (LBL), and 
;                Joe Hennawi (UC Berkeley)
;-
;------------------------------------------------------------------------------
FUNCTION FIRE_WSTRUCT_ld, hdr, LINELIST=linelist, REID_FILE=REID_FILE

;;-----------
;;  Line list directory
calib_path =  GETENV('FIRE_DIR') + '/LowDispersion/' ;GETENV('LONGSLIT_DIR') + '/calib/linelists/'
line_path = GETENV('FIRE_DIR') + '/LowDispersion/'
;;----------
;; Initial guess for central wavelength and dispersion
telescope = strtrim(sxpar(hdr, 'TELESCOP'))
instrument = strtrim(sxpar(hdr, 'INSTRUME'))
detector = strtrim(sxpar(hdr, 'DETECTOR'))
binning = long(strsplit(sxpar(hdr, 'BINNING'), ',', /extract))
mask    = strtrim(sxpar(hdr, 'SLITNAME'))
mswave = double(strtrim(sxpar(hdr, 'MSWAVE')))
telid     = strtrim(sxpar(hdr, 'TELID'), 2)


wstruct = create_struct('INSTRUMENT', instrument $
                        , 'REID', 0L             $
                        , 'AUTOID', 0L           $
                        , 'CRUDEID', 0L            $
                        , 'FWEIGHT', 0L            $
                        , 'PKWDTH', 0.0          $
                        , 'TOLER', 0.0           $
                        , 'NFIND', 0L            $
                        , 'DISP_GUESS', 0.0D     $
                        , 'FUNC', ' '            $
                        , 'NORD_CRUDE', 0L       $
                        , 'SIGREJ_CRUDE', 0.0    $ 
                        , 'PSIG_CRUDE', 0.0      $ 
                        , 'PSIG', fltarr(6)      $
                        , 'MXOFF', fltarr(6)     $  
                        , 'SIGREJ', fltarr(6)    $
                        , 'NORD', fltarr(6)      $ 
                        , 'FLG_QUAL', lonarr(6)      $ 
                        , 'PSIG_REID', fltarr(6)     $
                        , 'MXOFF_REID', fltarr(6)    $  
                        , 'SIGREJ_REID', fltarr(6)   $
                        , 'NORD_REID', lonarr(6)     $
                        , 'FLG_QUAL_REID', lonarr(6) $ 
                        , 'LINELIST_AP', ' '     $
                        , 'LINELIST',    ' '     $
                        , 'NPANIC', 0L           $
                        , 'NORD_PANIC', 0L       $
                        , 'DR_WAVE', [0.0, 0.0]  $
                        , 'MXSHIFT', 0.0  $
                        , 'REID_FILE', ' ' $
                        , 'BIN_RATIO',1. $  ;; JXP -- For LRISb 
                        , 'RADIUS', 0.0 $
                        , 'SIG_WPIX', 0.0)

wstruct.radius = 3.0
wstruct.SIG_WPIX = 3.0
if (instrument EQ 'LRISBLUE') then begin
    mask = strtrim(sxpar(hdr, 'SLITNAME'))
    grism = strtrim(sxpar(hdr, 'GRISNAME'))
    CASE grism OF
        '1200/3400': begin
            ;; parameters for arc peak finding
            wstruct.pkwdth        = 12.0/binning[1]
            wstruct.TOLER         = 4.0D/binning[1]
            wstruct.psig[0:4]     = [25.0, 18.0, 10.0, 7.0, 4.0]
            wstruct.MXOFF[0:4]    = [8.0, 8.0, 6.0, 6.0, 4.0]/binning[1]
            ;; paramters for wavelength soln fits
            wstruct.sigrej[0:4]   = [3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:4]     = [4L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:4] = [1L, 1L, 1L, 1L, 2L]
            wstruct.LINELIST      = line_path+'/lris_blue_1200.lst'
            wstruct.npanic        = 5L
            wstruct.nord_panic    = 2L
            wstruct.REID          = 1 ; reidentify using archived solutions
            wstruct.REID_FILE     = calib_path + '/lris_blue_1200.sav'
            wstruct.BIN_RATIO     = float(binning[1])/2L ;arxiv soln binned spectrally
        END
        '600/4000': begin
            ;; parameters for arc peak finding
            wstruct.pkwdth        = 12.0D/binning[1]
            wstruct.TOLER         = 4.0D/binning[1]
            wstruct.psig[0:4]     = [12.0, 8.0, 7.0, 4.0, 4.0]
            wstruct.mxoff[0:4]    = [7.0, 7.0, 5.0, 5.0, 5.0]/binning[1]
            ;; paramters for wavelength soln fits
            wstruct.sigrej[0:4]   = [3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:4]     = [4L, 5L, 5L, 5L, 5L]
            wstruct.FLG_QUAL[0:4] = [1L, 1L, 1L, 2L, 2L]
            wstruct.LINELIST      = line_path+'/lris_blue_600.lst'
            wstruct.npanic        = 5L
            wstruct.nord_panic    = 2L
            wstruct.REID          = 1 ; reidentify using archived solutions
            wstruct.REID_FILE     = calib_path + '/lris_blue_600.sav'
            wstruct.BIN_RATIO     = float(binning[1])/2L ;arxiv soln binned spectrally
        END
        '400/3400': begin
            wstruct.pkwdth      = 12.00/binning[1]
            wstruct.TOLER       = 4.0D/binning[1]
            wstruct.psig[0:4]   = [12.0, 8.0, 7.0, 4.0, 4.0]
            wstruct.mxoff[0:4]  = [7.0, 7.0, 5.0, 5.0, 5.0]/binning[1]
            wstruct.sigrej[0:4] = [3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:4]        = [4L, 5L, 5L, 5L, 5L]
            wstruct.FLG_QUAL[0:4]    = [1L, 1L, 1L, 2L, 2L]
            wstruct.LINELIST    = line_path+'/lris_blue_400.lst'
            wstruct.npanic      = 5L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1
            wstruct.REID_FILE   = calib_path + '/lris_blue_400.sav'
;           Archived solutions are binned 1 in spectral
            wstruct.BIN_RATIO   = binning[1]/1L
        END
        '300/5000': begin
            wstruct.pkwdth      = 12.0/binning[1] ; maxwidth of arc lines
            wstruct.TOLER       = 4.0D/binning[1]
            wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
            wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
            wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 2.5, 2.0]
            wstruct.nord[0:5]        = [4L, 4L, 5L, 5L, 6L, 6L]
            wstruct.FLG_QUAL[0:5]    = [2L, 2L, 2L, 3L, 3L, 3L]
            wstruct.LINELIST    = line_path+'lris_blue_300.lst'
            wstruct.npanic      = 5L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1 
            wstruct.REID_FILE   = calib_path + 'lris_blue_300.sav'
            wstruct.BIN_RATIO   = binning[1]/1L
        END
        ELSE: message, 'ERROR: Unknown grism'
    ENDCASE
ENDIF ELSE IF (instrument EQ 'LRIS') THEN BEGIN
    grating = strtrim(sxpar(hdr, 'GRANAME'))
    CASE grating OF
        '300/5000': BEGIN
            wstruct.pkwdth       = 6.0/binning[1]
            wstruct.TOLER        = 1.5/binning[1]
            wstruct.nfind        = 50L
            wstruct.disp_guess   = 2.45d0*binning[1]
            wstruct.FUNC  = 'CHEBY'
            wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
            wstruct.MXOFF[0:5]    = [6.0, 6.0, 6.0, 5.0, 3.0, 2.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST        = line_path+'/lris_red_300.lst'
            wstruct.npanic      = 20L
            wstruct.nord_panic  = 2L
            wstruct.dr_wave     = [0.9, 1.1]*(2.45d0)*binning[1]
            wstruct.AUTOID      = 1
        END
        '400/8500': BEGIN
            wstruct.pkwdth       = 6.0/binning[1]
            wstruct.TOLER        = 1.5/binning[1]
            wstruct.psig[0:5]  = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
            wstruct.MXOFF[0:5] = [6.0, 6.0, 6.0, 5.0, 3.0, 2.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST      = line_path+'/lris_red_300.lst'
            wstruct.npanic        = 10L
            wstruct.nord_panic    = 2L
            wstruct.REID        = 1
            wstruct.REID_FILE   = calib_path + $
              '/lris_red_400_grangle_24.32_multi.sav'
;;           These are AUTOID paramaters that should be used for a 
;;           a different grangle. 
            wstruct.nfind        = 30L
            wstruct.disp_guess   = 1.86d0*binning[1]
            wstruct.FUNC  = 'CHEBY'
            wstruct.dr_wave       = [0.9, 1.1]*(1.86d0)*binning[1]
            wstruct.AUTOID      = 1
        END
        '600/5000': BEGIN
            wstruct.pkwdth       = 6.0/binning[1]
            wstruct.TOLER        = 1.5/binning[1]
            wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
            wstruct.MXOFF[0:5]    = [5.0, 5.0, 5.0, 4.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:5]   = [4.0, 4.0, 4.0, 4.0, 4.0, 4.0]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 5L, 5L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST        = line_path+'/lris_red_600.lst'
            wstruct.npanic      = 20L
            wstruct.nord_panic  = 2L
            ;; Currently using autoid but we should archive multislit soln's
            wstruct.REID        = 0  
;            wstruct.REID_FILE   = calib_path + '/lris_red_600_7500.sav'
            ;; AUTOID paramaters
            wstruct.nfind        = 30L
            wstruct.disp_guess   = 1.28d0*binning[1]
            wstruct.FUNC  = 'CHEBY'
            wstruct.dr_wave     = [0.9, 1.1]*(1.28d0)*binning[1]
            wstruct.AUTOID      = 1
        END
        '600/7500': BEGIN
            wstruct.pkwdth       = 6.0/binning[1]
            wstruct.TOLER        = 1.5/binning[1]
            wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
            wstruct.MXOFF[0:5]    = [5.0, 5.0, 5.0, 4.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:5]   = [4.0, 4.0, 4.0, 4.0, 4.0, 4.0]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 5L, 5L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST        = line_path+'/lris_red_600.lst'
            wstruct.npanic      = 20L
            wstruct.nord_panic  = 2L
            ;; Currently using autoid but we should archive multislit soln's
            wstruct.REID        = 0  
            wstruct.REID_FILE   = calib_path + '/lris_red_600_7500.sav'
            ;; AUTOID paramaters
            wstruct.nfind        = 30L
            wstruct.disp_guess   = 1.28d0*binning[1]
            wstruct.FUNC  = 'CHEBY'
            wstruct.dr_wave     = [0.9, 1.1]*(1.28d0)*binning[1]
            wstruct.AUTOID      = 1
        END
        '600/10000': BEGIN
            wstruct.pkwdth       = 7.5/binning[1]
            wstruct.TOLER        = 2.0/binning[1]
            wstruct.FUNC  = 'CHEBY'
;            wstruct.psig[0] = 5.0
;            wstruct.MXOFF[0] = 3.0/binning[1]
;            wstruct.sigrej[0] = 3.0
;            wstruct.nord[0] = 5
;            wstruct.FLG_QUAL[0] = 2L
             wstruct.psig[0:5]     = [10.0, 10.0, 7.0, 7.0, 5.0, 5.0]
             wstruct.MXOFF[0:5]    = [3.0, 3.0, 3.0, 3.0, 2.0, 2.0]/binning[1]
             wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 2.5, 2.5]
             wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 5L, 5L]
             wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST        = line_path+'/lris_red_600.lst'
            wstruct.npanic      = 10L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 0
            wstruct.REID_FILE   = calib_path + $
              '/lris_red_600-10000_d680_mswave_8051.sav'
            wstruct.dr_wave     = [0.9, 1.1]*(1.28d0)*binning[1]
            wstruct.AUTOID      = 1
            wstruct.nfind        = 30L
            wstruct.disp_guess   = 1.28d0*binning[1]
        END
        '900/5500': BEGIN
            wstruct.pkwdth       = 6.0/binning[1]
            wstruct.TOLER        = 1.5/binning[1]
            wstruct.nfind        = 30L
            wstruct.disp_guess   = 0.85d0*binning[1]
            wstruct.FUNC  = 'CHEBY'
            wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
            wstruct.MXOFF[0:5]    = [6.0, 6.0, 6.0, 5.0, 3.0, 2.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST        = line_path+'/lris_red_600.lst'
            wstruct.REID        = 1
            wstruct.REID_FILE   = calib_path + '/lris_red_900_5500.sav'
            wstruct.npanic      = 20L
            wstruct.nord_panic  = 2L
            wstruct.dr_wave     = [0.9, 1.1]*(0.85d0)*binning[1]
            wstruct.AUTOID      = 1
        END
        '1200/7500': BEGIN
            wstruct.pkwdth        = 8.0/binning[1] ; width of peaks to center on
            wstruct.TOLER         = 3.0/binning[1] ; minimum centroiding diff
            wstruct.psig[0:5]     = [10.0, 7.0, 7.0, 5.0, 4.0, 3.0]
            wstruct.MXOFF[0:5]    = [4.0, 4.0, 4.0, 4.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [4L, 4L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:5] = [2L, 2L, 2L, 2L, 2L, 2L]
            wstruct.LINELIST      = line_path+'/lris_red_1200.lst'
            wstruct.npanic      = 4L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1
            wstruct.REID_FILE   = calib_path + $
              '/lris_red_1200-7500_d460_mswave_5137.sav'
            ;IF strmatch(mask, '*long*') AND mswave GT 6000.0D THEN BEGIN
;            wstruct.nfind       = 50L ; nlines for arc_pairs
;            wstruct.AUTOID      = 1L
;            wstruct.disp_guess   = 0.636d0*binning[1] ; dispersion
;            wstruct.dr_wave      = [0.95, 1.05]*(0.636d0)*binning[1]
;            wstruct.FUNC  = 'CHEBY' ; function to use for firts
            ;;ENDIF ELSE BEGIN
            ;;ENDELSE
        END
        ELSE: message, 'ERROR: Unknown grating'
    ENDCASE
ENDIF ELSE IF strmatch(instrument,'DEIMOS*') THEN BEGIN
    grating = strtrim(sxpar(hdr, 'GRATENAM'),2)
    gpos = strtrim(sxpar(hdr, 'GRATEPOS'),2)
    wave = sxpar(hdr,'G'+strtrim(gpos,2)+'TLTWAV')
    CASE grating OF
        '1200G': BEGIN
            wstruct.pkwdth       = 5.0/binning[1]
            wstruct.TOLER        = 2.0/binning[1]
            wstruct.nfind        = 50L
            wstruct.FUNC  = 'CHEBY'
            wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
            wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
            wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
            wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.npanic      = 15L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1
            wstruct.MXSHIFT     = 50
            if wave LT 7000 then begin
                wstruct.LINELIST  = line_path+'/deimos_blue.lst' 
                wstruct.REID_FILE   = calib_path + '/deimos_1200_blue.sav'
            endif else stop
        END
        ELSE: message, 'ERROR: Unknown grating'
    ENDCASE
endif else if (stregex(instrument,'.*kast.*',/boolean,/fold_case) eq 1) OR $
  (stregex(sxpar(hdr,'VERSION'),'kast*', /bool, /fold_case) EQ 1) then begin
    binning = [1,1]
    ;; Blue side
    if strtrim(sxpar(hdr,'SPSIDE'),2) EQ 'blue' OR $
      strmid(sxpar(hdr,'VERSION'),0,5) EQ 'kastb' then begin
        grating = strtrim(sxpar(hdr,'GRISM_N'),2)
        CASE grating OF
            '830/3460': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 2L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 50
                wstruct.LINELIST  = line_path+'/kast_blue3.lst' 
                wstruct.REID_FILE   = calib_path + '/kast_800_3460.sav'
            END
            '600/4310': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [6.0, 6.0, 6.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 4L, 5L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 2L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 200
                wstruct.LINELIST  = line_path+'/kast_blue2.lst' 
                wstruct.REID_FILE   = calib_path + '/kast_600_4310.sav'
            END
            '452/3306': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 2L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 50
                mjd_newccd = 2454534L 
                dateobs = sxpar(hdr, 'DATE-OBS')
                mjd = x_setjdate(strmid(dateobs,0,10))

                wstruct.LINELIST  = line_path+'/kast_blue3.lst' 
                if mjd LT 2454534L then $
                  wstruct.REID_FILE   = calib_path + '/kast_old_452_3306.sav'$
                else wstruct.REID_FILE   = calib_path + '/kast_452_3306.sav'
            END
            ELSE: message, 'ERROR: Unknown grating'
        ENDCASE
    endif else begin
        grating = strtrim(sxpar(hdr,'GRATNG_N'),2)
        CASE grating OF
            '1200/5000': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.disp_guess   = 1.17d0*binning[1]
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 3L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 50
                wstruct.LINELIST  = line_path+'/kast_red.lst' 
                wstruct.REID_FILE   = calib_path + '/kast_1200_5000.sav'
            END
            '1200 5000': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.disp_guess   = 1.17d0*binning[1]
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 3L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 50
                wstruct.LINELIST  = line_path+'/kast_red.lst' 
                wstruct.REID_FILE   = calib_path + '/kast_1200_5000.sav'
            END
            '600/5000': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.nfind        = 50L
                wstruct.disp_guess   = 2.35d0*binning[1]
                wstruct.FUNC  = 'CHEBY'
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 3L
;                wstruct.REID        = 1
;                wstruct.MXSHIFT     = 50
                wstruct.LINELIST  = line_path+'/kast_red.lst' 
                wstruct.AUTOID      = 1
                wstruct.dr_wave     = [0.9, 1.1]*(2.35d0)*binning[1]
;                wstruct.REID_FILE   = calib_path + '/kast_1200_5000.sav'
            END
            '600/7500': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.nfind        = 50L
                wstruct.disp_guess   = 2.32d0*binning[1]
                wstruct.FUNC  = 'CHEBY'
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 3L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 50
                wstruct.LINELIST  = line_path+'/kast_red.lst' 
;                wstruct.AUTOID      = 1
                wstruct.dr_wave     = [0.9, 1.1]*(2.32d0)*binning[1]
                wstruct.REID_FILE   = calib_path + '/kast_600_7500.sav'
            END
            ELSE: message, 'ERROR: Unknown grating'
        ENDCASE
    endelse
ENDIF ELSE IF strmatch(instrument,'DIS') THEN BEGIN
    binning = [1, 1]
    ;; Blue side
    grating = strtrim(sxpar(hdr, 'GRATING'), 2)
    CASE detector OF
        'blue': BEGIN
            wstruct.pkwdth        = 3.0/binning[1]
            wstruct.TOLER         = 2.0/binning[1]
            wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 6.0, 5.0]
            wstruct.MXOFF[0:5]    = [10.0, 8.0, 6.0, 4.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 2.5]
            wstruct.nord[0:5]     = [4L, 4L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.npanic        = 10L
            wstruct.nord_panic    = 2L
            wstruct.REID          = 1
            wstruct.MXSHIFT       = 200
            wstruct.LINELIST      = line_path+'/dis_henear.lst'
            wstruct.REID_FILE     = calib_path + '/dis_B400.sav'
        END
        'red': BEGIN
            wstruct.pkwdth        = 5.0/binning[1]
            wstruct.TOLER         = 2.0/binning[1]
            wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
            wstruct.MXOFF[0:5]    = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 5L, 5L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.npanic        = 10L
            wstruct.nord_panic    = 2L
            wstruct.REID          = 1
            wstruct.MXSHIFT       = 200
            wstruct.LINELIST  = line_path+'/dis_henear.lst'
            wstruct.REID_FILE   = calib_path + '/dis_R300.sav'
        END
        ELSE: message, 'ERROR: Unknown grating'
     ENDCASE
 ENDIF ELSE IF strmatch(telescope, 'kp4m') THEN BEGIN
    binning = [1, 1]
    grating = strcompress(sxpar(hdr, 'DISPERSE'), /rem)
    CASE grating OF
       'KPC10A': BEGIN
          wstruct.pkwdth        = 3.0/binning[1]
          wstruct.TOLER         = 2.0/binning[1]
          wstruct.psig[0:1]     = [10.0, 10.0]
          wstruct.MXOFF[0:1]    = [2.0, 2.0]/binning[1]
          wstruct.sigrej[0:1]   = [2.0, 2.0]
          wstruct.nord[0:1]     = [5L, 5L]
          wstruct.FLG_QUAL[0:1] = [2L, 2L]
          wstruct.FUNC  = 'CHEBY'
          wstruct.npanic        = 10L
          wstruct.nord_panic    = 2L
          wstruct.REID          = 1
          wstruct.MXSHIFT       = 200
          wstruct.LINELIST      = line_path+'/dis_henear.lst'
          wstruct.REID_FILE     = calib_path + '/kp4m_KPC10A.sav'
        END
        ELSE: message, 'ERROR: Unknown grating'
     ENDCASE
 ENDIF ELSE IF strmatch(telid, '200') THEN BEGIN
    binning = [1, 1]
    fpa = strtrim(sxpar(hdr, 'FPA'), 2) 
    CASE fpa OF
       'DBSP_BLUE': BEGIN
          wstruct.pkwdth        = 5.0/binning[1]
          wstruct.TOLER         = 3.0/binning[1]
          wstruct.psig[0:5]     = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
          wstruct.MXOFF[0:5]    = [4.0, 4.0, 3.0, 3.0, 3.0, 3.0]/binning[1]
          wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
          wstruct.nord[0:5]     = [5L, 5L, 5L, 5L, 5L, 5L]
          wstruct.FLG_QUAL[0:5] = [3L, 3L, 3L, 3L, 3L, 3L]
          wstruct.npanic        = 10L
          wstruct.nord_panic    = 2L
          wstruct.REID          = 1
          wstruct.MXSHIFT       = 200
          wstruct.LINELIST      = line_path+'/p200_DBSP_blue.lst'
          wstruct.REID_FILE     = calib_path + '/p200_DBSP_blue.sav'
        END
        'DBSP_RED': BEGIN
           wstruct.pkwdth       = 5.0/binning[1]
            wstruct.TOLER        = 3.0/binning[1]
            wstruct.nfind        = 40L
            wstruct.disp_guess   = 2.12d0*binning[1]
            wstruct.FUNC  = 'CHEBY'
            wstruct.psig[0:5]     = [10.0, 7.0, 7.0, 5.0, 5.0, 5.0]
            wstruct.MXOFF[0:5]    = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST        = line_path+'/dis_henear.lst'
            wstruct.npanic      = 10L
            wstruct.nord_panic  = 2L
            wstruct.dr_wave     = [0.9, 1.1]*(2.46d0)*binning[1]
            wstruct.AUTOID      = 1
        END
        ELSE: message, 'ERROR: Unknown grating'
    ENDCASE
ENDIF ELSE IF strmatch(telescope, 'mmt') THEN BEGIN
    grating = strtrim(sxpar(hdr, 'DISPERSE'))
    detector = sxpar(hdr, 'DETECTOR')
    n1 = strtrim(sxpar(hdr, 'NAXIS1'))
    n2 = strtrim(sxpar(hdr, 'NAXIS2'))
    binning = [round(2700./n1), 1]
    ;; MMT red channel?
    IF strmatch(detector, '*ccd34*') THEN grating = 'R150'
    case grating of 
        '800GPM': begin
            wstruct.pkwdth      = 8.00/binning[1]
            wstruct.TOLER       = 4.0D/binning[1]
            wstruct.psig[0:4]   = [12.0, 8.0, 7.0, 4.0, 4.0]
            wstruct.mxoff[0:4]  = [7.0, 7.0, 5.0, 5.0, 5.0]/binning[1]
            wstruct.sigrej[0:4] = [3.0, 3.0, 2.5, 2.5, 2.5]
            wstruct.nord[0:4]        = [3L, 3L, 3L, 3L, 3L]
            wstruct.FLG_QUAL[0:4]    = [1L, 1L, 1L, 2L, 2L]
            wstruct.npanic      = 5L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1
            if strmatch(sxpar(hdr,'COMPLAMP'), 'HeNeAr*') then begin
                wstruct.REID_FILE   = calib_path + '/mmt_blue_800_nocu.sav'
                wstruct.LINELIST    = line_path+'/henear_low.lst'
            endif else begin
                wstruct.REID_FILE   = calib_path + '/mmt_blue_800.sav'
                wstruct.LINELIST    = line_path+'/cuhenear.lst'
            endelse
        END
        '300GPM': begin
            wstruct.pkwdth      = 8.0D/binning[1]
            wstruct.TOLER       = 3.0D/binning[1]
            wstruct.FUNC  = 'CHEBY'
            wstruct.psig[0:4]   = [12.0, 8.0, 7.0, 5.0, 5.0]
            wstruct.mxoff[0:4]  = [7.0, 7.0, 3.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:4] = [3.0, 3.0, 2.5, 2.5, 2.5]
            wstruct.nord[0:4]        = [4L, 4L, 5L, 5L, 5L]
            wstruct.FLG_QUAL[0:4]    = [1L, 1L, 1L, 2L, 2L]
            wstruct.npanic      = 10L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1
            wstruct.REID_FILE   = calib_path + '/mmt_blue_300.sav'
            wstruct.LINELIST    = line_path + '/mmt_blue_300.lst'
        END
        'R150': BEGIN
            wstruct.pkwdth        = 7.0/binning[1]
            wstruct.TOLER         = 2.0/binning[1]
            wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
            wstruct.MXOFF[0:5]    = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 5L, 5L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.npanic        = 10L
            wstruct.nord_panic    = 2L
            wstruct.REID          = 1
            wstruct.MXSHIFT       = 200
            wstruct.LINELIST  = line_path+'/dis_henear.lst'
            wstruct.REID_FILE   = calib_path + '/mmt_red_150.sav'
        END
        ELSE: STOP
    ENDCASE
ENDIF ELSE IF strmatch(instrument, '*IMACS*') THEN BEGIN
    grating = strtrim(sxpar(hdr, 'DISPERSR'))
    binning = long(strsplit(sxpar(hdr, 'BINNING'), 'x', /extract))
    CASE grating OF
        'Gri-200-15.0': BEGIN
            wstruct.pkwdth        = 7.0/binning[1]
            wstruct.TOLER         = 2.0/binning[1]
            wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
            wstruct.MXOFF[0:5]    = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 5L, 5L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.npanic        = 10L
            wstruct.nord_panic    = 2L
            wstruct.REID          = 1
            wstruct.MXSHIFT       = 200
            wstruct.LINELIST  = line_path+'/dis_henear.lst'
            wstruct.REID_FILE   = calib_path + '/imacs_red_grism200.sav'
        END
        ELSE: STOP
    ENDCASE
ENDIF ELSE IF strmatch(telescope, 'Gemini-North') THEN BEGIN
    grating = strtrim(sxpar(hdr, 'GRATING'))
    binning = [1, 1] ;; assumes ny=1024 is default binning (need to generalize)
    wave = float(strcompress(sxpar(hdr, 'GRWLEN'), /rem))
    CASE grating OF
        'R150+_G5306': BEGIN
            wstruct.pkwdth        = 11.0/binning[1] ;; allow some fat lines
            wstruct.TOLER         = 3.0/binning[1]
            wstruct.psig[0:2]     = [5.0, 5.0, 4.0]
            wstruct.MXOFF[0:2]    = [3.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:2]   = [3.0, 3.0, 3.0]
            wstruct.nord[0:2]     = [3L, 3L, 3L]
            wstruct.FLG_QUAL[0:2] = [1L, 2L, 2L] ;, 1L, 1L, 2L, 2L, 2L]
            wstruct.npanic        = 4L
            wstruct.nord_panic    = 2L
            wstruct.REID          = 1
            wstruct.MXSHIFT       = 200
            wstruct.LINELIST  = line_path+'/GMOS_CuAr_R150.lst'
            wstruct.REID_FILE   = calib_path + '/GMOS_CuAr_R150_multislit.sav'
            wstruct.radius      = 2.0 ;; special for GMOS microslits
        END
        ELSE: STOP
    ENDCASE

 ENDIF ELSE IF (instrument EQ 'FIRE') then begin
    wstruct.pkwdth        = 12.0
    wstruct.TOLER         = 4.0D
    wstruct.psig[0:4]     = [25.0, 18.0, 10.0, 7.0, 4.0]
    wstruct.MXOFF[0:4]    = [8.0, 8.0, 6.0, 6.0, 4.0]
    ;; paramters for wavelength soln fits
    wstruct.sigrej[0:4]   = [3.0, 3.0, 3.0, 3.0, 3.0]
    wstruct.nord[0:4]     = [4L, 4L, 4L, 4L, 4L]
    wstruct.FLG_QUAL[0:4] = [1L, 1L, 1L, 1L, 2L]
    wstruct.LINELIST  = line_path+'/NeNeAr_fire.lst'
    wstruct.npanic        = 5L
    wstruct.nord_panic    = 2L
    wstruct.REID          = 1   ; reidentify using archived solutions
    wstruct.mxshift       = 50.0
    wstruct.REID_FILE     = getenv("FIRE_DIR")+'/LowDispersion/near_fed6.sav'
    wstruct.func          = 'CHEBY'
    wstruct.nord           = 6
ENDIF ELSE BEGIN
    message, 'ERROR: Unknown instrument!'
ENDELSE

if keyword_set(LINELIST) then wstruct.LINELIST = line_path+LINELIST
if keyword_set(REID_FILE) then wstruct.REID_FILE = calib_path+REID_FILE

RETURN, wstruct
END

