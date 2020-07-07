;;;;;;;;;;;;;;;;;;
pro gridsed, inputcounts, inputcounterrs=inputcounterrs, sed=sed


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

nsteps=5

counts_array=make_array(6,/FLOAT,VALUE=!values.f_nan) 

sedflux_lambda=interpol(stellarsed,sedwave, lambda)

for f=0,5 do counts_array[f]=tsum(lambda,filter_array[f,*]*sedflux_lambda*lambda/hc)

sedarray=stellarsed
chisqarray=total(((inputcounts-counts_array)/inputcounterrs)^2.0)

plot, sedwave, stellarsed

; note q is not an integer so that q/nsteps won't be an integer

sedflux=transpose(pickedfactors)*inputcounts[0:5]
sedflux=[sedflux,sedflux[5]]
newflux=sedflux
range=[4.0, 4, 2, 2, 0.2, 0.2]


stellarfluxdensityfactors=[6.2, 8.50, 4.0, 1.63, 1.472, 2.614] *10.0^(-16)
starthighfluxdensityfactors=[6.59, 9.38, 4.52, 1.79, 1.61, 3.10] *10.0^(-16)*2.0
startlowfluxdensityfactors=[0.01, 2.90, 0.48, 1.18, 1.07, 2.18] *10.0^(-16)/2.0

range=starthighfluxdensityfactors-startlowfluxdensityfactors




for l=0.0, nsteps-1 do begin
	newflux[0]=[ ( startlowfluxdensityfactors[0] + l/nsteps*(range[0]) )*inputcounts[0] ]
	for m=0.0, nsteps-1 do begin
		newflux[1]=[ ( startlowfluxdensityfactors[1] + m/nsteps*(range[1]) )*inputcounts[1] ]
		for n=0.0, nsteps-1 do begin
			newflux[2]=[ ( startlowfluxdensityfactors[2] + n/nsteps*(range[2]) )*inputcounts[2] ]
			for o=0.0, nsteps-1 do begin
				newflux[3]=[ ( startlowfluxdensityfactors[3] + o/nsteps*(range[3]) )*inputcounts[3] ]
				for p=0.0, nsteps-1 do begin
					newflux[4]=[ ( startlowfluxdensityfactors[4] + p/nsteps*(range[4]) )*inputcounts[4] ]
					for q=0.0,nsteps-1 do begin
						newflux[5]=[ ( startlowfluxdensityfactors[5] + q/nsteps*(range[5]) )*inputcounts[5] ]
						newflux[6]=newflux[5]
						sedflux_lambda=interpol(newflux,sedwave, lambda)
						for ff=0,5 do counts_array[ff]=tsum(lambda,filter_array[ff,*]*sedflux_lambda*lambda/hc)
						chisqarray=[chisqarray,total(((inputcounts-counts_array)/inputcounterrs)^2.0)]
						sedarray=[[sedarray], [newflux]]
						oplot, sedwave, newflux
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

plothist, chisqarray, xrange=[0,5]


fraction=0.2

chisqdif = 4.0

while min(chisqarray) gt 0.5 do begin

low=1.0-fraction
high=1.0+fraction

lowflux=dblarr(7)
highflux=dblarr(7)
for f=0,6 do lowflux[f]=min( sedarray[f,where(chisqarray lt min(chisqarray + chisqdif))] )
for f=0,6 do highflux[f]=max( sedarray[f,where(chisqarray lt min(chisqarray + chisqdif))] )



for l=0.0, nsteps-1 do begin
	newflux[0]=lowflux[0]+((highflux[0]-lowflux[0])*l/nsteps)
	for m=0.0, nsteps-1 do begin
		newflux[1]=lowflux[1]+((highflux[1]-lowflux[1])*m/nsteps)
		for n=0.0, nsteps-1 do begin
			newflux[2]=lowflux[2]+((highflux[2]-lowflux[2])*n/nsteps)
			for o=0.0, nsteps-1 do begin
				newflux[3]=lowflux[3]+((highflux[3]-lowflux[3])*o/nsteps)
				for p=0.0, nsteps-1 do begin
					newflux[4]=lowflux[4]+((highflux[4]-lowflux[4])*p/nsteps)
					for q=0.0,nsteps-1 do begin
						newflux[5]=lowflux[5]+((highflux[5]-lowflux[5])*q/nsteps)
						newflux[6]=newflux[5]
						sedflux_lambda=interpol(newflux,sedwave, lambda)
						for ff=0,5 do counts_array[ff]=tsum(lambda,filter_array[ff,*]*sedflux_lambda*lambda/hc)
						chisqarray=[chisqarray,total(((inputcounts-counts_array)/inputcounterrs)^2.0)]
						sedarray=[[sedarray], [newflux]]
						oplot, sedwave, newflux
					endfor
				endfor
			endfor
		endfor
	endfor
endfor
	allbest=where(chisqarray eq min(chisqarray))
print, min(chisqarray)
plothist, chisqarray, xrange=[0,5]
;fraction=fraction*0.8

chisqdif=chisqdif/2.0

endwhile


print, 'final stop'
stop
end
