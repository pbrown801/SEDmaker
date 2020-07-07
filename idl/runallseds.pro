;
pro runstandardsed

restore,  filename='spectrumfilenames.sav'

;;;;;;;;;;;;;; standard sed

sedstandard_array=dblarr(nspectra,9)
allmagsstandard_array=fltarr(nspectra,8)
allcountsstandard_array=fltarr(nspectra,8)
intfluxstandard_array=fltarr(nspectra)
muvfluxstandard_array=fltarr(nspectra)
nuvfluxstandard_array=fltarr(nspectra)
optfluxstandard_array=fltarr(nspectra)

filtereffwavelength=[1930,2200,2600,3450,4350,5460]

	sedstandard_wave=[1600,filtereffwavelength[0:5],6000,8000]

nomfluxfactors=[7.024, 9.179, 4.027, 1.666, 1.391, 2.728]*10.0^(-16.0)
stellarfluxfactors=[7.024, 9.179, 4.027, 1.666, 1.391, 2.728]*10.0^(-16.0)
filterzpts=[17.35, 16.82, 17.49, 18.34, 19.11, 17.89]


for n=0,nspectra-1 do begin
	spectrum=spectrumfiles[n]

	print, n, ' of ', nspectra, ' ', spectrum

	fluxsed=allcounts_array[n,0:5]*stellarfluxfactors[0:5]	
	

	sedstandard_array[n,*]=[0,reform(fluxsed),0,0]


	pjb_uvotspec_all, [transpose(sedstandard_wave),(sedstandard_array[n,*])], all=all

	allmagsstandard_array[n,*]=all.mag_array
	allcountsstandard_array[n,*]=all.counts_array
	intfluxstandard_array[n]=all.intflux
	muvfluxstandard_array[n]=all.muvflux
	nuvfluxstandard_array[n]=all.nuvflux
	optfluxstandard_array[n]=all.optflux

endfor

	save, filename='standardbolomags.sav', allmagsstandard_array,allcountsstandard_array,intfluxstandard_array,muvfluxstandard_array,nuvfluxstandard_array,optfluxstandard_array,sedstandard_array, sedstandard_wave


end 

pro runaveragesed

restore,  filename='spectrumfilenames.sav'


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;; sed average


sedaverage_array=dblarr(nspectra,9)
allmagsaverage_array=fltarr(nspectra,8)
allcountsaverage_array=fltarr(nspectra,8)
intfluxaverage_array=fltarr(nspectra)
muvfluxaverage_array=fltarr(nspectra)
nuvfluxaverage_array=fltarr(nspectra)
optfluxaverage_array=fltarr(nspectra)

filtereffwavelength=[1930,2200,2600,3450,4350,5460]

	sedaverage_wave=[1600,filtereffwavelength[0:5],6000,8000]

nomfluxfactors=[7.024, 9.179, 4.027, 1.666, 1.391, 2.728]*10.0^(-16.0)
stellarfluxfactors=[7.024, 9.179, 4.027, 1.666, 1.391, 2.728]*10.0^(-16.0)
filterzpts=[17.35, 16.82, 17.49, 18.34, 19.11, 17.89]


for n=0,nspectra-1 do begin
	spectrum=spectrumfiles[n]

	print, n, ' of ', nspectra, ' ', spectrum

	fluxsed=allcounts_array[n,0:5]*stellarfluxfactors[0:5]	
	

	sedaverage_array[n,*]=[fluxsed[0],reform(fluxsed),fluxsed[5],fluxsed[5]]


	pjb_uvotspec_all, [transpose(sedaverage_wave),(sedaverage_array[n,*])], all=all

	allmagsaverage_array[n,*]=all.mag_array
	allcountsaverage_array[n,*]=all.counts_array
	intfluxaverage_array[n]=all.intflux
	muvfluxaverage_array[n]=all.muvflux
	nuvfluxaverage_array[n]=all.nuvflux
	optfluxaverage_array[n]=all.optflux

endfor

	save, filename='averagebolomags.sav', allmagsaverage_array,allcountsaverage_array,intfluxaverage_array,muvfluxaverage_array,nuvfluxaverage_array,optfluxaverage_array,sedaverage_array, sedaverage_wave



end 

pro runbbsed

restore,  filename='spectrumfilenames.sav'


;;;;;;;;;;;;;;;   blackbody 



sedbb_array=dblarr(nspectra,9)
allmagsbb_array=fltarr(nspectra,8)
allcountsbb_array=fltarr(nspectra,8)
intfluxbb_array=fltarr(nspectra)
muvfluxbb_array=fltarr(nspectra)
nuvfluxbb_array=fltarr(nspectra)
optfluxbb_array=fltarr(nspectra)


for n=0,nspectra-1 do begin
	spectrum=spectrumfiles[n]

	print, n, ' of ', nspectra, ' ', spectrum

uvotcounts=allcounts_array[n,0:5]

fitbbcounts, uvotcounts=uvotcounts, besttemp=besttemp, lambda=lambda, planckflux_lambda=planckflux_lambda

	pjb_uvotspec_all, [transpose(lambda),transpose(planckflux_lambda)], all=all

	allmagsbb_array[n,*]=all.mag_array
	allcountsbb_array[n,*]=all.counts_array
	intfluxbb_array[n]=all.intflux
	muvfluxbb_array[n]=all.muvflux
	nuvfluxbb_array[n]=all.nuvflux
	optfluxbb_array[n]=all.optflux
	sedbb_wave=lambda
endfor

	save, filename='bbbolomags.sav', allmagsbb_array,allcountsbb_array,intfluxbb_array,muvfluxbb_array,nuvfluxbb_array,optfluxbb_array,sedbb_array, sedbb_wave



end 

pro rungridsed

restore,  filename='spectrumfilenames.sav'


;;;;;;;;;;;;;;; best fit sed


sedgrid_array=dblarr(nspectra,7)
fluxdensityfactorssedgrid_array=fltarr(nspectra,8)
allmagssedgrid_array=fltarr(nspectra,8)
allcountssedgrid_array=fltarr(nspectra,8)
intfluxsedgrid_array=fltarr(nspectra)
muvfluxsedgrid_array=fltarr(nspectra)
nuvfluxsedgrid_array=fltarr(nspectra)
optfluxsedgrid_array=fltarr(nspectra)

for n=0,nspectra-1 do begin
	spectrum=spectrumfiles[n]

	print, n, ' of ', nspectra, ' ', spectrum

uvotcounts=allcounts_array[n,0:5]


	gridsed, uvotcounts, bestsed=bestsed, sedwave=sedwave

	sedgrid_array[n,*]=bestsed
	sedgridwavelengths=sedwave
	sedgrid_wave=sedwave

	pjb_uvotspec_all, [transpose(sedgrid_wave),transpose(bestsed)], all=all

	allmagssedgrid_array[n,*]=all.mag_array
	allcountssedgrid_array[n,*]=all.counts_array
	intfluxsedgrid_array[n]=all.intflux
	muvfluxsedgrid_array[n]=all.muvflux
	nuvfluxsedgrid_array[n]=all.nuvflux
	optfluxsedgrid_array[n]=all.optflux

endfor


	save, filename='sedgridbolomags.sav', allmagssedgrid_array,allcountssedgrid_array,intfluxsedgrid_array,muvfluxsedgrid_array,nuvfluxsedgrid_array,optfluxsedgrid_array,sedgrid_array, sedgrid_wave, sedgridwavelengths



end 

pro runpickspec

restore,  filename='spectrumfilenames.sav'


;;;;;;;;;;; mangled spectrum


sedps_array=dblarr(nspectra,641)
allmagsps_array=fltarr(nspectra,8)
allcountsps_array=fltarr(nspectra,8)
intfluxps_array=fltarr(nspectra)
muvfluxps_array=fltarr(nspectra)
nuvfluxps_array=fltarr(nspectra)
optfluxps_array=fltarr(nspectra)
specpick_array=strarr(nspectra)

for n=0,nspectra-1 do begin
	spectrum=spectrumfiles[n]

	print, 'mangling spectra for ', n, ' of ', nspectra, ' ', spectrum

colordif=((allmags_array[n,1]-allmags_array[n,3])  - (allmags_array[*,1]-allmags_array[*,3] ))^2.0 + ((allmags_array[n,3]-allmags_array[n,5]) -  (allmags_array[*,3]-allmags_array[*,5] ))^2.0

	rank=sort(colordif)

specpick_array[n]=spectrumfiles[rank[4]]


uvotcounts=allcounts_array[n,0:5]
templatespectrum=specpick_array[n]


warpsed, templatespectrum, warpedspectrum, uvotcounts=uvotcounts

	pjb_uvotspec_all, [transpose(lambda),transpose(warpedspectrum)], all=all

	allmagsps_array[n,*]=all.mag_array
	allcountsps_array[n,*]=all.counts_array
	intfluxps_array[n]=all.intflux
	muvfluxps_array[n]=all.muvflux
	nuvfluxps_array[n]=all.nuvflux
	optfluxps_array[n]=all.optflux
;	sedps_array[n,*]=warpedspectrum[1,*]

endfor

	save, filename='pickspecmags.sav', allmagsps_array,allcountsps_array,intfluxps_array,muvfluxps_array,nuvfluxps_array,optfluxps_array,sedps_array, specpick_array


end

pro runallseds




end 

pro runallseds



runstandardsed
runaveragesed
runbbsed
rungridsed
runpickspec


print, 'final stop'
stop
end
