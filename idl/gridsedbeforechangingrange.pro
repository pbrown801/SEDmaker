;;;;;;;;;;;;;;;;;;
pro gridsed, inputcounts, sed, inputcounterrs=inputcounterrs


;  10,000 K black body
;   inputcounts=[5.00085e+06,  4.11938e+06,  8.76916e+06,  2.29198e+07,  1.97022e+07,  7.37948e+06]

; vega
; inputcounts=[8.97155e+06 , 5.51979e+06 , 9.52015e+06 , 2.18502e+07 , 4.41138e+07 , 1.43316e+07 ]
; 92a
; inputcounts=[4.54193   ,   1.81368  ,17.0897, 152.233 , 322.275 , 127.479 ]



IF KEYWORD_SET(inputcounterrs) THEN errorflag=1 else inputcounterrs=inputcounts*0.1

if n_elements(inputcounterrs) eq 1 then inputcounterrs=[inputcounterrs,inputcounterrs,inputcounterrs,inputcounterrs,inputcounterrs,inputcounterrs]

; SN Ia  inputmags=[15.714011 , 16.179403 , 14.405791 , 12.884552 ,12.839177 , 12.626133 ]

zeropoints=[17.35, 16.82, 17.49, 18.34, 19.11, 17.89]
;inputmagcounts=10.0^((zeropoints-inputmags)/2.5)

h = 1D*6.62606957E-27
c = 1D*2.99792458E18 ; in angstroms per second
hc = 1D*1.986E-8
hc = h*c
e = 1D*2.71828
k = 1D*1.38E-16

;read in filter curves
filtercurves=["$SNSCRIPTS/UVW2_2010.txt","$SNSCRIPTS/UVM2_2010.txt","$SNSCRIPTS/UVW1_2010.txt","$SNSCRIPTS/U_UVOT.txt","$SNSCRIPTS/B_UVOT.txt","$SNSCRIPTS/V_UVOT.txt"]
readcol,"$SNSCRIPTS/V_UVOT.txt",lambda,V_EA,/silent
filter_array=fltarr(n_elements(filtercurves),n_elements(lambda) )

for f=0,n_elements(filtercurves)-1 do begin
	readcol,filtercurves[f],filter_wave,filter_ea,/silent
	filter_array[f,*]=filter_ea
endfor


restore,  filename='spectrumfilenames.sav'

colordif=fltarr(nspectra)
for i=0,nspectra-1 do colordif[i]=total(abs((inputcounts/inputcounts[5]-allcounts_array/allcounts_array[5]))/(inputcounterrs/inputcounts[5]))


rank=sort(colordif)

picked=rank[0]
pickedfactors=FLUXDENSITYFACTORS_ARRAY[picked,0:5]

; default flux conversion factors
stellarfluxdensityfactors=[6.2, 8.50, 4.0, 1.63, 1.472, 2.614] *10.0^(-16)
vegaeffwave=[2030,2230,2590,3500,4330,5400,1950,2500]


sedwave=[1600,vegaeffwave[1:4],6000,8000]

stellarsed = [inputcounts*stellarfluxdensityfactors[0:5],inputcounts[5]*stellarfluxdensityfactors[5]]

nsteps=11

counts_array=make_array(6,/FLOAT,VALUE=!values.f_nan) 

sedflux_lambda=interpol(stellarsed,sedwave, lambda)

for f=0,5 do counts_array[f]=tsum(lambda,filter_array[f,*]*sedflux_lambda*lambda/hc)

sedarray=stellarsed
chisqarray=total(abs(inputcounts-counts_array)/inputcounterrs)

plot, sedwave, stellarsed

; note q is not an integer so that q/nsteps won't be an integer

sedflux=transpose(pickedfactors)*inputcounts[0:5]
sedflux=[sedflux,sedflux[5]]
newflux=sedflux
range=[4.0, 4, 2, 2, 0.2, 0.2]

;;; modify this range so nothing goes below zero
print, modify the range


;;; use minimum flux conversion





;range=starthighfluxdensityfactors-startlowfluxdensityfactors

;	sedflux= [ ( startlowfluxdensityfactors + q/nsteps*(range) )*inputcounts[0:5],( startlowfluxdensityfactors[5] + q/nsteps*(range[5]) )*inputcounts[5] ]






stop

for l=0.0, nsteps-1 do begin
	newflux[0]=sedflux[0]+(l-5)/nsteps*sedflux[0]*range[0]
	for m=0.0, nsteps-1 do begin
		newflux[1]=sedflux[1]+(m-5)/nsteps*sedflux[1]*range[1]
		for n=0.0, nsteps-1 do begin
			newflux[2]=sedflux[2]+(n-5)/nsteps*sedflux[2]*range[2]
			for o=0.0, 2 do begin
				newflux[3]=sedflux[3]+(o-1)/nsteps*sedflux[3]*range[3]
				for p=0.0, 2 do begin
					newflux[4]=sedflux[4]+(p-1)/nsteps*sedflux[4]*range[4]
					;print, newflux[4]
					for q=0.0,2 do begin

						newflux[5]=sedflux[5]+(q-1)/nsteps*sedflux[5]*range[5]
						newflux[6]=newflux[5]
						;print, q, newflux

						sedflux_lambda=interpol(newflux,sedwave, lambda)

						for ff=0,5 do counts_array[ff]=tsum(lambda,filter_array[ff,*]*sedflux_lambda*lambda/hc)

						chisqarray=[chisqarray,total(abs(inputcounts-counts_array)/inputcounterrs)]
						sedarray=[[sedarray], [newflux]]
						oplot, sedwave, newflux
						;print, newflux[5]
					endfor
				endfor
			endfor
		endfor
	endfor
endfor
	allbest=where(chisqarray eq min(chisqarray))

print, min(chisqarray)
oplot, sedwave, sedarray[*,allbest], thick=2

;plot, sedarray[1,*], chisqarray
;;;;;;;


print, 'final stop'
stop
end
