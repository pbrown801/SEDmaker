pro makebolosedplots5, spectrum

figurename='sed5plots_'+spectrum+'.eps'


	pjb_uvotspec_all, spectrum, mag_array=mag_array, counts_array=counts_array, lambda=lambda, filter_array=filter_array, sp_wave=sp_wave, sp_flux=sp_flux, smooth_wave=smooth_wave, smooth_flux=smooth_flux

ylog=alog10(max(smooth_flux))
top=ceil(max(smooth_flux)/10.0^floor(ylog))+1

ymax=top

nplots=5
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
nxticks=22
nyticks=top

x1 = margin*8.8/xsize
x2 = x1 + a*8.8/xsize
xc = x2 + wall*8.8/xsize
y1 = margin*8.8/ysize
y2 = y1 + b*8.8/ysize
xrange=[1600,6000]
yrange=[0,ymax]
yrange1=[0,ymax]
yrange2=[0,ymax]


SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=12, bits_per_pixel=8, /color

x=0
plot, sp_wave,sp_flux/10.0^floor(ylog), /noerase, $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle='Wavelength [Angstroms]',  $
; ytitle='Flux Density ', charsize=1.0,  $
 ytitle='Flux Density / 10^'+strtrim(floor(ylog),1), charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=[' ',' ','2000', ' ', ' ', ' ', ' ','3000', ' ', ' ', ' ', ' ','4000', ' ', ' ', ' ', ' ','5000', ' ', ' ', ' ', ' ','6000']

;;; pick the 10th closest spectrum

restore,  filename='spectrumfilenames.sav'

colordif=((mag_array[1]-mag_array[3])  - (allmags_array[*,1]-allmags_array[*,3] ))^2.0 + ((mag_array[3]-mag_array[5]) -  (allmags_array[*,3]-allmags_array[*,5] ))^2.0

	rank=sort(colordif)

templatespectrum=spectrumfiles[rank[9]]

uvotcounts=counts_array

warpsed, uvotcounts=uvotcounts, templatespectrum, warpedspectrum

oplot, lambda, warpedspectrum/10.0^floor(ylog), thick=2

vend=interpol(warpedspectrum,lambda,6000.0)


xpoly=[1600.0, lambda[0:440] ,6000.0,6000.0,1600.0]
ypoly=[0.0, reform(warpedspectrum[0:440]), vend, 0.0 ,0.0]/10.0^floor(ylog)
cgcolorfill, xpoly, ypoly, color='grey'

cgoplot, sp_wave,sp_flux/10.0^floor(ylog), linestyle=0

xyouts, 4400, ymax-1, 'Warped Spectrum'

plot, sp_wave,sp_flux/10.0^floor(ylog), /noerase, $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle='Wavelength [Angstroms]',  $
; ytitle='Flux Density ', charsize=1.0,  $
 ytitle='Flux Density / 10^'+strtrim(floor(ylog),1), charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=[' ',' ','2000', ' ', ' ', ' ', ' ','3000', ' ', ' ', ' ', ' ','4000', ' ', ' ', ' ', ' ','5000', ' ', ' ', ' ', ' ','6000']

x=1
plot, sp_wave,sp_flux/10.0^floor(ylog), /noerase, $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=' ',   ytitle='Flux Density / 10^'+strtrim(floor(ylog),1), charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

gridsed, uvotcounts, bestsed=bestsed, sedarray=sedarray, chisqarray=chisqarray, sedcount_array=sedcount_array, sigmaone=sigmaone, sedwave=sedwave


oplot, sedwave[0:5], bestsed[0:5]/10.0^floor(ylog)

xpoly=[sedwave[0:5] , 6000.0,1600.0]
ypoly=[reform(bestsed[0:5]), 0.0 ,0.0 ]/10.0^floor(ylog)
cgcolorfill, xpoly, ypoly, color='grey'

oplot, sp_wave,sp_flux/10.0^floor(ylog), linestyle=0, thick=2

cgoplot, sedwave[0:5], bestsed[0:5]/10.0^floor(ylog), psym=-14


xyouts, 4400, ymax-1, 'Best-fit SED'


x=2
plot, sp_wave,sp_flux/10.0^floor(ylog), /noerase, $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=' ',   ytitle='Flux Density / 10^'+strtrim(floor(ylog),1), charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

uvotcounts=counts_array
fitbbcounts, uvotcounts=uvotcounts, counterrs=counterrs, besttemp=besttemp, lambda=lambda, planckflux_lambda=planckflux_lambda, plancklum=plancklum

;stop

oplot, lambda, planckflux_lambda/10.0^floor(ylog)

xpoly=[lambda[0:440] , 6000.0, 1600.0]
ypoly=[reform(planckflux_lambda[0:440]), 0.0, 0.0 ]/10.0^floor(ylog)
cgcolorfill, xpoly, ypoly, color='grey'

oplot, sp_wave,sp_flux/10.0^floor(ylog), linestyle=0, thick=2


cgoplot, lambda, planckflux_lambda/10.0^floor(ylog)

xyouts, 4400, ymax-1, 'Blackbody'

x=3
plot, sp_wave,sp_flux/10.0^floor(ylog), /noerase, $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=xtitle6,   ytitle='Flux Density / 10^'+strtrim(floor(ylog),1), charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

filtereffwavelength=[1930,2200,2600,3450,4350,5460]

sedstandard_wave=[1600,filtereffwavelength[0:5],6000,8000]
sedaverage_wave=[1600,filtereffwavelength[0:5],6000,8000]
nomfluxfactors=[7.024, 9.179, 4.027, 1.666, 1.391, 2.728]*10.0^(-16.0)
stellarfluxfactors=[7.024, 9.179, 4.027, 1.666, 1.391, 2.728]*10.0^(-16.0)
filterzpts=[17.35, 16.82, 17.49, 18.34, 19.11, 17.89]

	fluxsed=counts_array[0:5]*stellarfluxfactors[0:5]	
	sedaverage_flux=[fluxsed[0],reform(fluxsed),fluxsed[5],fluxsed[5]]
	
	sedstandard_flux=[0,reform(fluxsed),0,0]

oplot, sedaverage_wave, sedaverage_flux/10.0^floor(ylog)


xpoly=[sedaverage_wave[0:7] ,6000.0, 1600.0]
ypoly=[reform(sedaverage_flux[0:7]) ,0.0, 0.0 ]/10.0^floor(ylog)
cgcolorfill, xpoly, ypoly, color='grey'


oplot, sp_wave,sp_flux/10.0^floor(ylog), linestyle=0, thick=2

cgoplot, sedaverage_wave, sedaverage_flux/10.0^floor(ylog), psym=-14

xyouts, 4400, ymax-1, 'Standard-Flat'

x=4
plot, sp_wave,sp_flux/10.0^floor(ylog), /noerase, $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=xtitle6,   ytitle='Flux Density / 10^'+strtrim(floor(ylog),1), charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

oplot, sedstandard_wave, sedstandard_flux/10.0^floor(ylog)


xpoly=[sedstandard_wave[0:6] , 6000.0, 1600.0]
ypoly=[reform(sedstandard_flux[0:6]) ,0.0, 0.0 ]/10.0^floor(ylog)
cgcolorfill, xpoly, ypoly, color='grey'


oplot, sp_wave,sp_flux/10.0^floor(ylog), linestyle=0, thick=2

cgoplot, sedstandard_wave, sedstandard_flux/10.0^floor(ylog), psym=-14

xyouts, 4400, ymax-1, 'Standard-Zero'

device, /close
SET_PLOT, 'X'
$open sed5plots_vega.flm.eps

print, 'final stop'
stop
end

