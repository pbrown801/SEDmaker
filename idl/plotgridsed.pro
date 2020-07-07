;;;;;;;;;;;;;;;;;;
pro plotgridsed, inputspectrum

;inputspectrum='sn1992A-19920124-iue-ctio.flm'
;inputspectrum='vega.flm'

pjb_uvotspec_all, inputspectrum, all=all

inputcounts=all.counts_array



IF KEYWORD_SET(inputcounterrs) THEN errorflag=1 else inputcounterrs=inputcounts*0.1

if n_elements(inputcounterrs) eq 1 then inputcounterrs=[inputcounterrs,inputcounterrs,inputcounterrs,inputcounterrs,inputcounterrs,inputcounterrs]

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


;restore,  filename='spectrumfilenames.sav'
;colordif=fltarr(nspectra)
;for i=0,nspectra-1 do colordif[i]=total(abs((inputcounts/inputcounts[5]-allcounts_array/allcounts_array[5]))/(inputcounterrs/inputcounts[5]))
;rank=sort(colordif)
;picked=rank[0]
;pickedfactors=FLUXDENSITYFACTORS_ARRAY[picked,0:5]
;sedflux=transpose(pickedfactors)*inputcounts[0:5]
;sedflux=[sedflux,sedflux[5]]
;newflux=sedflux

; default flux conversion factors
stellarfluxdensityfactors=[6.2, 8.50, 4.0, 1.63, 1.472, 2.614] *10.0^(-16)
vegaeffwave=[2030,2230,2590,3500,4330,5400,1950,2500]


sedwave=[1600,vegaeffwave[1:4],6000,8000]

stellarsed = [inputcounts*stellarfluxdensityfactors[0:5],inputcounts[5]*stellarfluxdensityfactors[5]]

nsteps=4

counts_array=make_array(6,/FLOAT,VALUE=!values.f_nan) 

sedflux_lambda=interpol(stellarsed,sedwave, lambda)

for f=0,5 do counts_array[f]=tsum(lambda,filter_array[f,*]*sedflux_lambda*lambda/hc)

sedarray=stellarsed
newflux=stellarsed
chisqarray=total(((inputcounts-counts_array)/inputcounterrs)^2.0)



nplots=1
; from http://www.iluvatar.org/~dwijn/idlfigures
!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
; the default size is given in centimeters
; 8.8 is made to match a journal column width
xsize = 8.8
wall = 0.03
margin=0.12
a = xsize/8.8 - (margin + wall)
b = a * 2d / (1 + sqrt(5))

ysize = (margin + nplots*(b + wall ) )*xsize
ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a

x1 = margin*8.8/xsize
x2 = x1 + a*8.8/xsize
xc = x2 + wall*8.8/xsize
y1 = margin*8.8/ysize
y2 = y1 + b*8.8/ysize

ytitle='Flux Density'


figurename=inputspectrum+'gridsed.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=12, bits_per_pixel=8, /color

cgplot, all.sp_wave, alog10(all.sp_flux), color='red', thick=3, position=[x1,y1,x2,y2], $
xtitle='Wavelength [Angstroms]',   ytitle=ytitle, charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
ystyle=1, xrange=[1600,6000], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues


oplot, sedwave, alog10(stellarsed)



stellarfluxdensityfactors=[6.2, 8.50, 4.0, 1.63, 1.472, 2.614] *10.0^(-16)
starthighfluxdensityfactors=[6.59, 9.38, 4.52, 1.79, 1.61, 3.10] *10.0^(-16)*1.0
startlowfluxdensityfactors=[0.01, 2.90, 0.48, 1.18, 1.07, 2.18] *10.0^(-16)/1.0


iteration=0

while min(chisqarray) gt 0.8 and iteration lt 5 do begin


chisqdif = max([min(chisqarray)*2.0,1.0])


lowflux=dblarr(7)
highflux=dblarr(7)

if iteration eq 0 then for f=0,5 do lowflux[f]  =startlowfluxdensityfactors[f]*inputcounts[f]
if iteration eq 0 then for f=0,5 do highflux[f]=starthighfluxdensityfactors[f]*inputcounts[f]


if iteration gt 0 then for f=0,6 do lowflux[f]=min( sedarray[f,where(chisqarray lt min(chisqarray + chisqdif))] )*0.8
if iteration gt 0 then for f=0,6 do highflux[f]=max( sedarray[f,where(chisqarray lt min(chisqarray + chisqdif))] )*1.2



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
						oplot, sedwave, alog10(newflux)
					endfor
				endfor
			endfor
		endfor
	endfor
endfor

	allbest=where(chisqarray eq min(chisqarray))

print, min(chisqarray)
cgoplot, sedwave, alog10(sedarray[*,allbest]), color='blue', thick=2
cgoplot, all.sp_wave, alog10(all.sp_flux), color='red', thick=3
;plot, sedarray[1,*], chisqarray
;;;;;;;

;plothist, chisqarray, xrange=[0,5]


chisqdif=chisqdif/2.0

iteration=iteration+1.0

endwhile


bestsed=sedarray[*,allbest]


device, /close
SET_PLOT, 'X'
$open *gridsed.eps 

stop

end
