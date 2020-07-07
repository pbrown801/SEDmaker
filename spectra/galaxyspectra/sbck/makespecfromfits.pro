pro makespecfromfits, fitsspectrum


print, fitsspectrum
;fitsspectrum='ugca410.fits'

flux=mrdfits(fitsspectrum,0,header)

;  from header
;CRVAL1  =              2792.46 
;CDELT1  =                 5.33 

CRVAL1=sxpar(header, 'CRVAL1')
CDELT1=sxpar(header, 'CDELT1')

wave=fltarr(n_elements(flux))
for q=0,n_elements(flux)-1 do wave[q]=CRVAL1+q*CDELT1
plot, wave, flux, xrange=[1600,6000], xstyle=1


openw,lun,fitsspectrum+'.dat',/get_lun

for q=0,n_elements(flux)-1 do printf, lun, wave[q], flux[q]

free_lun,lun


end