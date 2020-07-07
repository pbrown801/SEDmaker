;;;;;;;;;;;;;;;;;;
pro makesedcounts, inputcounts, sed, inputcounterrs=inputcounterrs sedarray=sedarray


;  10,000 K black body
;   inputcounts=[5.00085e+06,  4.11938e+06,  8.76916e+06,  2.29198e+07,  1.97022e+07,  7.37948e+06]

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
starthighfluxdensityfactors=[6.59, 9.38, 4.52, 1.79, 1.61, 3.10] *10.0^(-16)*1.2
startlowfluxdensityfactors=[0.01, 2.90, 0.48, 1.18, 1.07, 2.18] *10.0^(-16)/1.2

sed6flux=inputcounts*stellarfluxdensityfactors[0:5]

filtereffwavelengths=[1930,2200,2600,3450,4350,5460]
vegaeffwave=[2030,2230,2590,3500,4330,5400,1950,2500]

inputflux=inputcounts*stellarfluxdensityfactors[0:5]


sedwave=[1600,vegaeffwave[1:4],6000,8000]

stellarsed = [inputcounts*stellarfluxdensityfactors[0:5],inputcounts[5]*stellarfluxdensityfactors[5]]

nsteps=10

counts_array=make_array(6,/FLOAT,VALUE=!values.f_nan) 
shortchisqarray=make_array(nsteps,/FLOAT,VALUE=!values.f_nan) 
shortsedarray=make_array(7,nsteps,/FLOAT,VALUE=!values.f_nan) 
;  an array for all of the seds and chi square values will be made on the fly and added to
range=(starthighfluxdensityfactors-startlowfluxdensityfactors)*inputcounts


sedflux_lambda=interpol(stellarsed,sedwave, lambda)

for f=0,5 do counts_array[f]=tsum(lambda,filter_array[f,*]*sedflux_lambda*lambda/hc)

sedarray=stellarsed
chisqarray=total(abs(inputcounts-counts_array)/inputcounterrs)

plot, sedwave, stellarsed

; note q is not an integer so that q/nsteps won't be an integer

;;; now start with the stellar seds and vary each filter one at a time in this order
filterorder=[1,5,4,3,2,1,0,3,4,5,4,3,2,1,0,1,5,4,3,2,1,0,5,4,3,2,1,0,5,4,3,2,1,0,5,4,3,2,1,0,5,4,3,2,1,0,5,4,3,2,1,0]


;	sedflux=stellarfluxdensityfactors*inputcounts[0:5]
	sedflux=transpose(pickedfactors)*inputcounts[0:5]
	sedflux=[sedflux,sedflux[5]]

for f=0, n_elements(filterorder)-1 do begin
filter=filterorder[f]


	if f eq 16 then range=range/5.0
	if f eq 22 then range=range/5.0
	if f eq 32 then range=range*5.0

	for q=0.0,nsteps-1 do begin

		filterflux=sedflux[filter] + (q/nsteps*(range[filter]) )*inputcounts[filter]
		sedflux[filter]=filterflux
		; to keep the last point equal to v
		sedflux[6]=sedflux[5]
		if f eq 15 then sedflux[0]=filterflux
		if f eq 15 then sedflux[2]=filterflux
;	print, 'f', f, 'q', q

		sedflux_lambda=interpol(sedflux,sedwave, lambda)

		for ff=0,5 do counts_array[ff]=tsum(lambda,filter_array[ff,*]*sedflux_lambda*lambda/hc)

		shortchisqarray[q]=total(abs(inputcounts-counts_array)/inputcounterrs)
		shortsedarray[*,q]=sedflux
		chisqarray=[chisqarray,total(abs(inputcounts-counts_array)/inputcounterrs)]
		sedarray=[[sedarray], [sedflux]]
oplot, sedwave, sedflux
wait, 0.1

	endfor

;	if f eq 7 then sedflux=stellarfluxdensityfactors*inputcounts[0:5]
	if f eq 7 then sedflux=transpose(pickedfactors)*inputcounts[0:5]
	if f eq 7 then sedflux=[sedflux,sedflux[5]]

	best=where(shortchisqarray eq min(shortchisqarray))
	if f ge 16 then sedflux=shortsedarray[*,best]
	allbest=where(chisqarray eq min(chisqarray))
	allbest=allbest[n_elements(allbest)-1]
	if f ge 22 then sedflux=sedarray[*,allbest]
	if f ge 22 then range=sedarray[*,allbest]*0.2
	if f ge 28 then range=sedarray[*,allbest]*0.1
	if f ge 34 then range=sedarray[*,allbest]*0.05
	if f ge 40 then range=sedarray[*,allbest]*0.1



endfor

print, min(chisqarray)
;plot, sedarray[1,*], chisqarray
;;;;;;;


print, 'final stop'
stop
end
