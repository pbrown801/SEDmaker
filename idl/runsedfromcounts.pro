

pro runsed7

restore,  filename='spectrumfilenames.sav'


sed7_array=dblarr(nspectra,7)
fluxdensityfactors7_array=fltarr(nspectra,8)
allmags7_array=fltarr(nspectra,8)
allcounts7_array=fltarr(nspectra,8)
intflux7_array=fltarr(nspectra)
muvflux7_array=fltarr(nspectra)
nuvflux7_array=fltarr(nspectra)
optflux7_array=fltarr(nspectra)

for n=0,nspectra-1 do begin
	spectrum=spectrumfiles[n]

	print, n, ' of ', nspectra, ' ', spectrum

	makesed7, allmags_array[n,*],  sed7

	sed7_array[n,*]=sed7[1,*]
	sed7wavelengths=sed7[0,*]
	sed7_wave=sed7[0,*]

	pjb_uvotspec_all, sed7, mag_array=mag_array, counts_array=counts_array, lambda=lambda, filter_array=filter_array, countspec_array=countspec_array,  fluxdensityfactors=fluxdensityfactors, flatflux=flatflux,  intflux=intflux, muvflux=muvflux, nuvflux=nuvflux, optflux=optflux

	allmags7_array[n,*]=mag_array
	allcounts7_array[n,*]=counts_array
	intflux7_array[n]=intflux
	muvflux7_array[n]=muvflux
	nuvflux7_array[n]=nuvflux
	optflux7_array[n]=optflux

endfor


	save, filename='sed7bolomags.sav', allmags7_array,allcounts7_array,intflux7_array,muvflux7_array,nuvflux7_array,optflux7_array,sed7_array, sed7_wave, sed7wavelengths



print, 'last stop'
stop
end
