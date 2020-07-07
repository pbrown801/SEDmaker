;;;;;;;;;;;;;;;;;;
pro makesed7, inputmags, sed


; SN Ia  inputmags=[15.714011 , 16.179403 , 14.405791 , 12.884552 ,12.839177 , 12.626133 ]

zeropoints=[17.35, 16.82, 17.49, 18.34, 19.11, 17.89]
inputmagcounts=10.0^((zeropoints-inputmags)/2.5)

; default flux conversion factors
stellarfluxdensityfactors=[6.2, 8.50, 4.0, 1.63, 1.472, 2.614] *10.0^(-16)

sed6flux=inputmagcounts*stellarfluxdensityfactors[0:5]

filtereffwavelengths=[1930,2200,2600,3450,4350,5460]

inputflux=inputmagcounts*stellarfluxdensityfactors[0:5]

readcol,"$SNSCRIPTS/V_UVOT.txt",lambda,V_EA,/silent

;print, ' start at blue end with m2 fixing the flux '
;print, ' '

nsteps=20.0
f=1
range=0.4

diff=fltarr(nsteps)
for q=0.0,nsteps-1 do begin

	flux=(1.0-range)*sed6flux[f] + q/nsteps*(2.0*range)*sed6flux[f]
	newwavelengths=[1600,8000]
	newflux=[flux,flux]
	newfluxspectrum=interpol(newflux,newwavelengths, lambda)
	spectrum=[transpose(lambda),transpose(newfluxspectrum)]
	pjb_uvotspec_all, spectrum, mag_array=mag_array
	diff[q]=(mag_array[f]-inputmags[f])
endfor

best=where(abs(diff) eq min(abs(diff)))
bestflux=(1.0-range)*sed6flux[f] + float(best[0])/nsteps*(2.0*range)*sed6flux[f]
flux=bestflux

m2flux=bestflux

;;;;;;;;
f=4


diff=fltarr(nsteps)
for q=0.0,nsteps-1 do begin

	flux=(1.0-range)*sed6flux[f] + q/nsteps*(2.0*range)*sed6flux[f]
	newwavelengths=[1600,8000]
	newflux=[flux,flux]
	newfluxspectrum=interpol(newflux,newwavelengths, lambda)
	spectrum=[transpose(lambda),transpose(newfluxspectrum)]
	pjb_uvotspec_all, spectrum, mag_array=mag_array
	diff[q]=(mag_array[f]-inputmags[f])
endfor

best=where(abs(diff) eq min(abs(diff)))
bestflux=(1.0-range)*sed6flux[f] + float(best[0])/nsteps*(2.0*range)*sed6flux[f]

bbflux=bestflux

sed7wavelengths=[1600,2200,2600,3450,4350,6000,8000]

sed7fluxes=[m2flux,m2flux,sed6flux[2],sed6flux[3],bbflux,sed6flux[5],sed6flux[5]]

newsed7fluxes=sed7fluxes

range=0.6

for i=0,n_elements(sed7wavelengths)-2 do begin
	h=n_elements(sed7wavelengths)-2-i
	newsed7fluxes=sed7fluxes

	diff=fltarr(nsteps)
	for q=0.0,nsteps-1 do begin

		flux=(1.0-range)*sed7fluxes[h] + q/nsteps*(2.0*range)*sed7fluxes[h]
		newsed7fluxes[h]=flux
		newsed7fluxes[6]=newsed7fluxes[5]
		newfluxspectrum=interpol(newsed7fluxes,sed7wavelengths, lambda)
		spectrum=[transpose(lambda),transpose(newfluxspectrum)]
		pjb_uvotspec_all, spectrum, mag_array=mag_array
		diff[q]=(total(abs(mag_array[h]-inputmags[h])))
	endfor

	best=where( diff eq min( diff ))
	bestflux=(1.0-range)*sed7fluxes[h] + float(best[0])/nsteps*(2.0*range)*sed7fluxes[h]

	sed7fluxes[h]=bestflux
	sed7fluxes[6]=sed7fluxes[5]
	newsed7fluxes=sed7fluxes
	pjb_uvotspec_all, [transpose(sed7wavelengths),transpose(sed7fluxes)], mag_array=mag_array

;	print, mag_array[0:5]

endfor


for h=0,n_elements(sed7wavelengths)-2 do begin

range=0.1


	diff=fltarr(nsteps)
	for q=0.0,nsteps-1 do begin

		flux=(1.0-range)*sed7fluxes[h] + q/nsteps*(2.0*range)*sed7fluxes[h]
		newsed7fluxes[h]=flux
		newsed7fluxes[6]=newsed7fluxes[5]
		newfluxspectrum=interpol(newsed7fluxes,sed7wavelengths, lambda)
		spectrum=[transpose(lambda),transpose(newfluxspectrum)]
		pjb_uvotspec_all, spectrum, mag_array=mag_array
		diff[q]=(total(abs(mag_array[h]-inputmags[h])))
	endfor

	best=where( diff eq min( diff ))
	bestflux=(1.0-range)*sed7fluxes[h] + float(best[0])/nsteps*(2.0*range)*sed7fluxes[h]

	sed7fluxes[h]=bestflux
	sed7fluxes[6]=sed7fluxes[5]
	newsed7fluxes=sed7fluxes
	pjb_uvotspec_all, [transpose(sed7wavelengths),transpose(sed7fluxes)], mag_array=mag_array

;	print, mag_array[0:5]


endfor


	newfluxspectrum=interpol(sed7fluxes,sed7wavelengths, lambda)
	spectrum=[transpose(lambda),transpose(newfluxspectrum)]


	pjb_uvotspec_all, spectrum, mag_array=mag_array, counts_array=counts_array, lambda=lambda, filter_array=filter_array, countspec_array=countspec_array,  fluxdensityfactors=fluxdensityfactors, flatflux=flatflux,  intflux=intflux, muvflux=muvflux, nuvflux=nuvflux, optflux=optflux




;print, inputmags
;print, mag_array[0:5]

sed=[transpose(sed7wavelengths),transpose(sed7fluxes)]



end


