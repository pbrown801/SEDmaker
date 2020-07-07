pro makekilonovased


;;;; create blackbody spectrum
h = 1D*6.626E-27
c = 1D*2.998E18
hc = 1D*1.986E-8
e = 1D*2.71828
k = 1D*1.38E-16
readcol,"$SNSCRIPTS/V_UVOT.txt",lambda,V_EA,/silent


uvotcounts=[0.147 , 0.094  , 0.576   , 2.814, 3.19,  1.64]


uvotcounts_err=[  0.035 ,  0.019 ,  0.069 ,  0.258,   0.19  ,  0.08]

gridsed, uvotcounts, bestsed=bestsed, sedarray=sedarray, chisqarray=chisqarray, sedcount_array=sedcount_array, sigmaone=sigmaone, sedwave=sedwave

fitbbcounts, uvotcounts=uvotcounts, counterrs=counterrs, besttemp=besttemp, lambda=lambda, planckflux_lambda=planckflux_lambda, plancklum=plancklum

print, besttemp
plot, lambda, planckflux_lambda

oplot, sedwave, bestsed, psym=-4

stellarfluxdensityfactors=[6.2, 8.50, 4.0, 1.63, 1.472, 2.614] *10.0^(-16)

stellarflux=stellarfluxdensityfactors*uvotcounts
oplot, sedwave, stellarflux, psym=-6
T1=8864.0
bbspec1=1L*2*h*c^2/(lambda^5*(e^(h*c/(lambda*k*T1))-1))*10^20.0

cgoplot, lambda, bbspec1*max(stellarflux)/max(bbspec1), color='blue'

spectrum=[transpose(lambda),transpose(bbspec1)]
pjb_uvotspec_all, spectrum, all=all

counts_array1=all.counts_array
speceffwavelength1=all.speceffwave
fluxspeceffwave1=all.fluxspeceffwave
factorspeceffwave1=all.factorspeceffwave
fluxdensityfactors1=all.factorvegaeffwave
vegaeffwavelength=all.vegaeffwave


cgoplot, vegaeffwavelength, fluxdensityfactors1*uvotcounts, psym=4, color='red', symsize=4

cgoplot, vegaeffwavelength, stellarfluxdensityfactors*uvotcounts, psym=4, color='dark green', symsize=5

T2=6000
bbspec2=1L*2*h*c^2/(lambda^5*(e^(h*c/(lambda*k*T1))-1))*10^20.0

cgoplot, lambda, bbspec2*uvotcounts[5]/counts_array2[5], color='blue'

spectrum=[transpose(lambda),transpose(bbspec2)]
pjb_uvotspec_all, spectrum, all=all

counts_array2=all.counts_array
speceffwavelength2=all.speceffwave
fluxspeceffwave2=all.fluxspeceffwave
factorspeceffwave2=all.factorspeceffwave
fluxdensityfactors2=all.factorvegaeffwave
vegaeffwavelength=all.vegaeffwave


cgoplot, vegaeffwavelength, fluxdensityfactors2*uvotcounts, psym=4, color='purple', symsize=4

cgoplot, lambda, bbspec1*uvotcounts[5]/counts_array1[5], color='purple'

stop

didn't go beyond here

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
xrange=[-15,35]
yrangetop=[0.0,0.55]
yrangebottom=[0.0,0.12]

xdata=[1,2,3,4]
ydata=[2,3,4,5]

nxticks=10
nyticks=8
ytitle='log (Flux Density)'

xrange=[1000,7000]
yrange=[(-19.0),(-12.0)]
yrange=[0,8]

figurename='comparekilonovamodels.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=12, bits_per_pixel=8, /color

plot, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='Wavelength [Angstroms]',   ytitle=ytitle, charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
yrange=[-1,8], ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks+1, ytickv=ytickvalues




cgoplot, lambda,alog10(bbspec2red), linestyle=2, color='blue'
cgoplot, lambda,alog10(bbspec1), color='red'

cgoplot, vegaeffwavelength[0:5], alog10(counts_array1[0:5]*fluxdensityfactors1), psym=5, symsize=1, symcolor='red'
cgoplot, vegaeffwavelength[0:5], alog10(counts_array1[0:5]*fluxdensityfactors2), psym=6, symsize=1, symcolor='blue'

xyouts, vegaeffwavelength[0]-100.0, alog10(counts_array1[0]*fluxdensityfactors2[0])-1.5, 'w2'
xyouts, vegaeffwavelength[1]-100.0, alog10(counts_array1[1]*fluxdensityfactors2[1])-1.2, 'm2'
xyouts, vegaeffwavelength[2]-100.0, alog10(counts_array1[2]*fluxdensityfactors2[2])-1.0, 'w1'
xyouts, vegaeffwavelength[3]-100.0, alog10(counts_array1[3]*fluxdensityfactors2[3])-1.0, 'u'
xyouts, vegaeffwavelength[4]-100.0, alog10(counts_array1[4]*fluxdensityfactors2[4])-1.0, 'b'
xyouts, vegaeffwavelength[5]-100.0, alog10(counts_array1[5]*fluxdensityfactors2[5])-1.0, 'v'


;al_legend, ['35000 K red phot', '35000 K red','2000 K phot','2000 K'],  psym=[6,0,5,0], $
;linestyle=[0,2,0,0], color=['blue','blue','red','red'], symsize=[1,0,1,0], $
;background_color='white', pos=[0.43,0.5], /norm, charsize=1.0, box=0
al_legend, ['35000 K red', '2000 K'],  psym=[0,0], linestyle=[2,0], color=['blue', 'red'], symsize=[0,0], background_color='white',$
pos=[0.43,0.54], /norm, charsize=1.0, box=0
al_legend, ['counts/sec x ff[35000 K,ex]', 'counts/sec x ff[2000 K]'],  psym=[6,5], $
color=['blue', 'red'], symsize=[1,1], $
background_color='white', pos=[0.43,0.4], /norm, charsize=1.0, box=0


device, /close
SET_PLOT, 'X'
$open comparebothmodelsfixed.eps 





print, 'final stop'
stop
end