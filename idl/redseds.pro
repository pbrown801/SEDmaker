pro redseds


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  3 PLOT  ;;;;;;;;;;;;;;;;;;;;
;

nplots=3
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
margin=0.2
a = xsize/8.8 - (margin + wall)
b = a * 2d / (1 + sqrt(5))

ysize = (margin + nplots*(b + wall ) )*xsize

; change to fit page size
xsize=18.0
ysize=23.0
a = xsize/8.8 - (1.5*margin + wall)
b=((ysize/8.8 - margin)/nplots - wall)/2
ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a

x1 = 1.5*margin*8.8/xsize
x2 = x1 + a*8.8/xsize
xc = x2 + wall*8.8/xsize
y1 = margin*8.8/ysize
y2 = y1 + b*8.8/ysize

xdata=[1,2,3,4]
ydata=[2,3,4,5]

figurename='redseds.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=12, bits_per_pixel=8, /color

q=0

;;;;;;;;;;;;;;;;;;;;;;;;;
xrange=[4000,6000]
nxticks=4

pjb_uvotspec_all, 'k4i.dat', counts_array=counts_array, sp_wave=sp_wave, sp_flux=sp_flux

plot, sp_wave, sp_flux,  /noerase, $
position=[x1,y1+(3.0-1)*b*xsize/ysize,x2,y1+(3.0+0)*b*xsize/ysize], $
xtitle=' ',   ytitle='Flux Density', charsize=1.5, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
$ yrange=[40,45.5],yticks=11, 
ystyle=1, xrange=xrange, xstyle=1,  $
xticks=nxticks, xtickv=xtickvalues, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1)

vegaeffwavelength=[2030,2231,2634,3501,4329,5402,2030,2634]
nomfluxfactors=[7.024, 9.179, 4.027, 1.666, 1.391, 2.728]*10.0^(-16.0)

oplot, vegaeffwavelength[3:5], nomfluxfactors[3:5]*counts_array[3:5], psym=4, symsize=2

;;;;;;;;;;;;;;;;;;;
pjb_uvotspec_all, '13784.dat' , counts_array=counts_array, sp_wave=sp_wave, sp_flux=sp_flux

plot, sp_wave, sp_flux, /noerase, $
position=[x1,y1+(2-1)*b*xsize/ysize,x2,y1+(2+0)*b*xsize/ysize], $
xtitle=' ',   ytitle='Flux Density', charsize=1.5, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
$ yrange=[5,50.0], /ylog, yticks=5, ystyle=1, 
xrange=xrange, xstyle=1, yrange=[0,2.0*10.0^(-15.0)],$
xticks=nxticks, xtickv=xtickvalues, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1)


oplot, vegaeffwavelength[3:5], nomfluxfactors[3:5]*counts_array[3:5], psym=4, symsize=2

;;;;;;;;;;;
restore, 'lbgspec_pjb.sav'
redshiftindex=375
pjb_uvotspec_all, [transpose(lc),transpose(reform(zlbg[redshiftindex,*])/10.0^20.0)], counts_array=counts_array, sp_wave=sp_wave, sp_flux=sp_flux


plot, sp_wave, sp_flux, /noerase, position=[x1,y1,x2,y1+(1+0)*b*xsize/ysize], $
xtitle='Observed Wavelength [Angstroms]',   ytitle='Flux Density', charsize=1.5, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
$ yrange=[0.1,100.0], /ylog, yticks=3,  ystyle=1, 
xrange=xrange, xstyle=1, yrange=[0.0,2.0*10.0^9.0], $
xticks=nxticks, xtickv=xtickvalues,  ytickv=ytickvalues

oplot, vegaeffwavelength[3:5], nomfluxfactors[3:5]*counts_array[3:5], psym=4, symsize=2

device, /close
SET_PLOT, 'X'
$open redseds.eps 







redshiftindex=365
pjb_uvotspec_all, [transpose(lc),transpose(reform(zlbg[redshiftindex,*])/10.0^20.0)], counts_array=counts_array, sp_wave=sp_wave, sp_flux=sp_flux


plot, sp_wave, sp_flux, /noerase,  $
xtitle='Observed Wavelength [Angstroms]',   ytitle='Flux Density', charsize=1.5, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
$ yrange=[0.1,100.0], /ylog, yticks=3,  ystyle=1, 
xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues,  ytickv=ytickvalues

oplot, vegaeffwavelength[3:5], nomfluxfactors[3:5]*counts_array[3:5], psym=4, symsize=2


;;;;;
redshiftindex=355
pjb_uvotspec_all, [transpose(lc),transpose(reform(zlbg[redshiftindex,*])/10.0^20.0)], counts_array=counts_array, sp_wave=sp_wave, sp_flux=sp_flux


plot, sp_wave, sp_flux, /noerase,  $
xtitle='Observed Wavelength [Angstroms]',   ytitle='Flux Density', charsize=1.5, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
$ yrange=[0.1,100.0], /ylog, yticks=3,  ystyle=1, 
xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues,  ytickv=ytickvalues

oplot, vegaeffwavelength[3:5], nomfluxfactors[3:5]*counts_array[3:5], psym=4, symsize=2


;;;;;
redshiftindex=365
pjb_uvotspec_all, [transpose(lc),transpose(reform(zlbg[redshiftindex,*])/10.0^20.0)], counts_array=counts_array, sp_wave=sp_wave, sp_flux=sp_flux


plot, sp_wave, sp_flux,  $
xtitle='Observed Wavelength [Angstroms]',   ytitle='Flux Density', charsize=1.5, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
$ yrange=[0.1,100.0], /ylog, yticks=3,  ystyle=1, 
xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues,  ytickv=ytickvalues

oplot, vegaeffwavelength[3:5], nomfluxfactors[3:5]*counts_array[3:5], psym=4, symsize=2


;;;;;




;;;;;
redshiftindex=375
pjb_uvotspec_all, [transpose(lc),transpose(reform(zlbg[redshiftindex,*])/10.0^20.0)], counts_array=counts_array, sp_wave=sp_wave, sp_flux=sp_flux


plot, sp_wave, sp_flux,  $
xtitle='Observed Wavelength [Angstroms]',   ytitle='Flux Density', charsize=1.5, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
$ yrange=[0.1,100.0], /ylog, yticks=3,  ystyle=1, 
xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues,  ytickv=ytickvalues

oplot, vegaeffwavelength[3:5], nomfluxfactors[3:5]*counts_array[3:5], psym=4, symsize=2


;;;;;




stop
end