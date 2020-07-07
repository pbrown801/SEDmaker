pro makeboloredsedplots, v
;make this so it reads in the spectra and plots from the spectrum array with the accompanying seds


;restore, 'vega_92a_bolocompare.sav'

;oldpath=!path
;olddir=!dir
restore, filename='vega_92a_redbolocompare_25.sav', verbose=0
;!path=oldpath
;!dir=olddir


spectrum_flux=allspectra[*,v]

spectrumvega='vega.dat'


readcol,spectrumfiles[v],sp_wave,sp_flux,/silent

filtereffwavelengths=[1930,2200,2600,3450,4350,5460]
nomfluxfactors=[7.024, 9.179, 4.027, 1.666, 1.391, 2.728]*10.0^(-16.0)
filterzpts=[17.35, 16.82, 17.49, 18.34, 19.11, 17.89]

sed13wavelengths=[1600,1930,2090,2200,2355,2600,3030,3450,3870,4350,4970,5460,8000]
sed6wavelengths=[1600,filtereffwavelengths,8000]

grbfluxdensityfactors=[6.2, 8.5, 4.0, 1.63, 1.472, 2.614] *10.0^(-16)
starfluxdensityfactors=[6.0, 7.5, 4.3, 1.5, 1.32, 2.61] *10.0^(-16)
filtereffwavelengths=[1930,2200,2600,3450,4350,5460]

;;;;;;;;;;;;;
;;Read in the filter Effective Area curves
readcol,"$SNSCRIPTS/V_UVOT.txt",lambda,V_EA,/silent
readcol,"$SNSCRIPTS/B_UVOT.txt",lambda,B_EA,/silent
readcol,"$SNSCRIPTS/U_UVOT.txt",lambda,U_EA,/silent
readcol,"$SNSCRIPTS/UVW1_2010.txt",lambda,W1_EA,/silent
readcol,"$SNSCRIPTS/UVM2_2010.txt",lambda,M2_EA,/silent
readcol,"$SNSCRIPTS/UVW2_2010.txt",lambda,W2_EA,/silent
readcol,"$SNSCRIPTS/UVW1-rc.txt",lambda,W1rc_EA,/silent
readcol,"$SNSCRIPTS/UVW2-rc.txt",lambda,W2rc_EA,/silent
;;;


nplots=2
; from http://www.iluvatar.org/~dwijn/idlfigures
!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
; the default size is given in centimeters
; 8.8 is made to match a journal column width
xsize = 8.8
wall = 0.05
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
xrange=[1500,6000]
nxticks=9

xticknames=['1600',' ','2000',' ',' ',' ',' ','3000',' ',' ',' ',' ','4000',' ',' ',' ',' ','5000',' ',' ',' ',' ','6000']

xticknames=['1500','2000','2500','3000','3500','4000','4500','5000','5500','6000']


specmax=max(spectrum_flux,/nan)

sedmax=max(allstellaravgendseddered[*,v],/nan)


x=string(max([specmax,sedmax],/nan))

exp=strsplit(x,'e', /extract)

yrange=[0,ceil(float(exp[0])) ]
yrange1=[0,4]
yrange2=[0,4]

nyticks=ceil(float(exp[0]))

ytitle='Flux Density'

yscale=10^float((-1.0)*exp[1]) ; 10.0^9.0
;  trim off peak value and set yaxis based on that

figurename=spectrumfiles[v]+'_deredsedplots.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=12, bits_per_pixel=8, /color

x=0

plot, lambda, spectrum_flux*yscale, /noerase, $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle='Wavelength [Angstroms]',   ytitle=ytitle, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=xticknames


xpoly=[sed7wavelengths[0:5] , 6000.0,1600.0]
ypoly=[all7seddered[0:5,v], 0.0,0.0 ]*yscale
cgcolorfill, xpoly, ypoly, color='grey'

oplot, lambda, spectrum_flux*yscale, linestyle=0, thick=2
oplot, sed7wavelengths, all7seddered[*,v]*yscale, psym=4

xyouts, 4800.0, ceil(float(exp[0]))-1, 'Floating SED'


x=1
plot, lambda, spectrum_flux*yscale, /noerase, $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=xtitle6,   ytitle=ytitle, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

;xpoly=[filterends[0:5] ,6000.0, 1600.0]
;ypoly=[allstellaravgendseddered[0:5,v] ,0.0, 0.0 ]*yscale
xpoly=[filterends[1:6] ,6000.0, 1600.0]
ypoly=[allstellaravgendseddered[1:6,v] ,0.0, 0.0 ]*yscale
cgcolorfill, xpoly, ypoly, color='grey'

oplot, lambda, spectrum_flux*yscale, linestyle=0, thick=2

oplot, filterends[1:6], allstellaravgendseddered[1:6,v]*yscale, psym=4


xyouts, 4800.0, ceil(float(exp[0]))-1, 'Standard'


x=0
plot, lambda, spectrum_flux*yscale, /noerase, $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle='Wavelength [Angstroms]',   ytitle=ytitle, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=xticknames

x=1
plot, lambda, spectrum_flux*yscale, /noerase, $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=xtitle6,   ytitle=ytitle, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)


device, /close
SET_PLOT, 'X'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



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
wall = 0.05
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

figurename=spectrumfiles[v]+'_deredsedsed7.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=12, bits_per_pixel=8, /color

x=0


plot, lambda, spectrum_flux*yscale, /noerase, $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle='Wavelength [Angstroms]',   ytitle=ytitle, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=xticknames


xpoly=[sed7wavelengths[0:5] , 6000.0,1600.0]
ypoly=[all7seddered[0:5,v], 0.0,0.0 ]*yscale
cgcolorfill, xpoly, ypoly, color='grey'

oplot, lambda, spectrum_flux*yscale, linestyle=0, thick=2
oplot, sed7wavelengths, all7seddered[*,v]*yscale, psym=4

xyouts, 4800.0, ceil(float(exp[0]))-1, 'Floating SED'



plot, lambda, spectrum_flux*yscale, /noerase, $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle='Wavelength [Angstroms]',   ytitle=ytitle, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=xticknames




device, /close
SET_PLOT, 'X'



figurename=spectrumfiles[v]+'_deredsedstandard.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=12, bits_per_pixel=8, /color

x=0

plot, lambda, spectrum_flux*yscale, /noerase, $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=xtitle6,   ytitle=ytitle, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=xticknames

;xpoly=[filterends[0:5] ,6000.0, 1600.0]
;ypoly=[allstellaravgendseddered[0:5,v] ,0.0, 0.0 ]*yscale
xpoly=[filterends[0:7] ,6000.0, 1600.0]
ypoly=[allstellaravgendseddered[0:7,v] ,0.0, 0.0 ]*yscale
cgcolorfill, xpoly, ypoly, color='grey'

oplot, lambda, spectrum_flux*yscale, linestyle=0, thick=2

oplot, filterends, allstellaravgendseddered[*,v]*yscale, psym=4


xyouts, 4800.0, ceil(float(exp[0]))-1, 'Standard'


plot, lambda, spectrum_flux*yscale, /noerase, $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=xtitle6,   ytitle=ytitle, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=xticknames



device, /close
SET_PLOT, 'X'




stop
end


