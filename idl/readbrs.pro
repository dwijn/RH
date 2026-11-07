PRO readBRS
;+
; NAME:
;	READBRS
;
; PURPOSE:
;	This routine reads the record structure for the background
;       opacities. For each wavelength index this structure specifies
;       the record number of the backgound opacities (all spatial locations)
;       for that wavelength.
;
; CATEGORY:
;	Data reduction
;
; CALLING SEQUENCE:
;       READBRS
;
; INPUTS:
;	brs.out (by default, hardwired).
;
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;	atmos.backgr_rn:  Record numbers of background opacity (and
;                         emissivity and scattering coefficient) for
;                         each wavelength index. The opacity data itself
;                         is usually stored in background.dat
;	atmos.backgrbytes:    Number of bytes per value in background.dat
;                         (4 = float, 8 = double). Newer versions of RH
;                         store the background opacities as floats and
;                         append this number as a final XDR int to
;                         brs.out. When it is absent (brs.out written by
;                         older versions) background.dat is in double
;                         precision and backgrbytes is set to 8.
;	atmos.backgretaexp:   Binary exponent of the scaling of the emissivities
;                         in background.dat: they are stored as
;                         eta * 2^backgretaexp (chi and scatt are not
;                         scaled). Follows backgrbytes as the last XDR int
;                         of brs.out; 0 when absent.
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
;       geometrycommon, atmoscommon
; SIDE EFFECTS:
; RESTRICTIONS:
;       The geometry of the atmosphere should be known (passed via COMMON block
;       geometrycommon), and the atmos structure should exist.
; PROCEDURE:
; EXAMPLE:
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Thu Jan  3 14:51:14 2008 --
;-

@geometry.common
@atmos.common

  openr, unit, "brs.out", /GET_LUN, /XDR

  atmosID = ''
  Nspace = 0L  &  Nspect = 0L
  readu, unit, atmosID, Nspace, Nspect

  hasline     = intarr(Nspect)
  ispolarized = intarr(Nspect)

  IF (atmos.moving OR tag_present(atmos, 'STOKES')) THEN $
   backgrrecno = lonarr(2*Nspect * geometry.Nrays) $
  ELSE $
   backgrrecno = lonarr(Nspect)

  readu, unit, hasline, ispolarized, backgrrecno

  ;; --- Optional trailing int: bytes per value in background.dat --- ;

  bgbytes = 8L  &  bgetaexp = 0L
  IF (NOT eof(unit)) THEN BEGIN
    bgbytes = 0L  &  bgetaexp = 0L
    readu, unit, bgbytes, bgetaexp
  ENDIF

  backgrflags = $
   replicate({hasline: hasline[0], ispolarized: ispolarized[0]}, Nspect)
  backgrflags.hasline     = hasline
  backgrflags.ispolarized = ispolarized

  ;; --- Add to atmos; update in place when readbrs is called again --- ;

  IF (tag_present(atmos, 'BACKGRFLAGS')) THEN BEGIN
    atmos.backgrflags = backgrflags
    atmos.backgrrecno = backgrrecno
    atmos.backgrbytes  = bgbytes
    atmos.backgretaexp = bgetaexp
  ENDIF ELSE $
   atmos = create_struct(atmos, 'backgrflags', backgrflags, $
                         'backgrrecno', backgrrecno, $
                         'backgrbytes', bgbytes, 'backgretaexp', bgetaexp)
  free_lun, unit
END
