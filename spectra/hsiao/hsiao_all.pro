pro specphot, spectrum, specmags, redfrac=redfrac

h = 1D*6.626E-27
c = 1D*2.998E18
hc = 1D*1.986E-8
e = 1D*2.71828
k = 1D*1.38E-16


;;Read in the filter Effective Area curves
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/V_UVOT.txt",lambda,V_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/B_UVOT.txt",lambda,B_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/U_UVOT.txt",lambda,U_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVW1.txt",lambda,W1_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVM2.txt",lambda,M2_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVW2.txt",lambda,W2_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVW1-rc.txt",lambda,W1rc_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVW2-rc.txt",lambda,W2rc_EA,/silent

;filters=fltarr(9,n_elements(lambda)
;filters[0,*]=lambda


;;Read in the Spectra and interpolate to the
;;Filter Curves

; check to see if it is a string (ie a filename) or an array
s=size(spectrum)
if (s[1] eq 7) then begin
	readcol,spectrum,sp_wave,sp_flux,/silent
endif else begin
	sp_wave=spectrum[0,*]
	sp_flux=spectrum[1,*]
endelse

flux=interpol(sp_flux,sp_wave,lambda)

v_counts=V_EA*flux*(10*lambda/hc)
b_counts=B_EA*flux*(10*lambda/hc)
u_counts=U_EA*flux*(10*lambda/hc)
w1_counts=W1_EA*flux*(10*lambda/hc)
m2_counts=M2_EA*flux*(10*lambda/hc)
w2_counts=W2_EA*flux*(10*lambda/hc)
w1rc_counts=W1rc_EA*flux*(10*lambda/hc)
w2rc_counts=W2rc_EA*flux*(10*lambda/hc)
;window, 0
;plot, lambda, m2_counts
;oplot, lambda, w2rc_counts

v_spec=-2.5*alog10(total(v_counts))+17.89
b_spec=-2.5*alog10(total(b_counts))+19.11
u_spec=-2.5*alog10(total(u_counts))+18.34
w1_spec=-2.5*alog10(total(w1_counts))+17.49
m2_spec=-2.5*alog10(total(m2_counts))+16.82
w2_spec=-2.5*alog10(total(w2_counts))+17.35
w1rc_spec=-2.5*alog10(total(w1rc_counts))+17.34
w2rc_spec=-2.5*alog10(total(w2rc_counts))+17.25

redfrac=[total(w1rc_counts)/total(w1_counts),total(w2rc_counts)/total(w2_counts)]
;specmags=[v_spec, b_spec, u_spec, w1_spec, m2_spec, w2_spec, w1rc_spec, w2rc_spec]

specmags=[w2_spec, m2_spec, w1_spec, u_spec, b_spec, v_spec, w1rc_spec, w2rc_spec]
;print, "specmags", specmags
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mangle, spectrum, snmags, z, mangledspectrum

h = 6.626E-27
c = 2.998E18
hc = 1.986E-8

w2_obs=snmags[0]
m2_obs=snmags[1]
w1_obs=snmags[2]
u_obs=snmags[3]  
b_obs=snmags[4]
v_obs=snmags[5]


;;Read in the filter Effective Area curves
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/V_UVOT.txt",lambda,V_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/B_UVOT.txt",lambda,B_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/U_UVOT.txt",lambda,U_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVW1.txt",lambda,W1_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVM2.txt",lambda,M2_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVW2.txt",lambda,W2_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVW1-rc.txt",lambda,W1rc_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVW2-rc.txt",lambda,W2rc_EA,/silent

; check to see if the spectrum is a string (ie a filename) or an array
s=size(spectrum)
if (s[1] eq 7) then begin
	readcol,spectrum,sp_wave,sp_flux,/silent
endif else begin
	sp_wave=spectrum[0,*]
	sp_flux=spectrum[1,*]
endelse

flux=interpol(sp_flux,sp_wave,lambda)
;plot, lambda, flux
sp_wave_z=sp_wave*(1+z)

flux_z=interpol(sp_flux,sp_wave_z,lambda)
flux_rest=interpol(sp_flux,sp_wave,lambda)

tempflux=sp_flux
x=0
	;;TYLER: Iterate to a convergence value, until each |ratio-1|<= 0.01
	filterratio=0
	while( (total(where(abs(filterratio-1) gt 0.01))  ) ne -1    and ( x lt 5) ) do begin
	x=x+1
	;;pass spectrum through effective area curves and convert to counts
	specphot, [(sp_wave_z), (tempflux)], specmags


	;;compare the synthetic magnitude to the observed magnitudes
	;;to create correction ratios
	v_ratio =2.512^(specmags[5]-v_obs)
	b_ratio =2.512^(specmags[4]-b_obs)
	u_ratio =2.512^(specmags[3]-u_obs)
	w1_ratio=2.512^(specmags[2]-w1_obs)
	m2_ratio=2.512^(specmags[1]-m2_obs)
	w2_ratio=2.512^(specmags[0]-w2_obs)


	;;create a correction term as a function of wavelength
	filtereffwavelength=[2231,2634,3501,4329,5402]
	filterratio=[m2_ratio, w1_ratio, u_ratio, b_ratio, v_ratio]
	filtereffwavelength=[1600, 2500,3000,3501,4329,5402, 8000]
	filterratio=[m2_ratio,m2_ratio, w1_ratio, u_ratio, b_ratio, v_ratio, v_ratio]
; this is from sne_mangle_new
	filtereffwavelength=[1600, 2231,2634,3501,4329,5402,8000]
	filterratio=[m2_ratio,m2_ratio, w1_ratio, u_ratio, b_ratio, v_ratio, v_ratio]
;print, filterratio
	responsecurve=interpol(filterratio,filtereffwavelength,sp_wave_z)
;plot, sp_wave_z, responsecurve
	tempflux=tempflux*responsecurve
	endwhile

	mangledspectrum=[(sp_wave_z),(tempflux)]
;window, 3
;plot, sp_wave_z, tempflux
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro kcor, restframespectrum, z, kcor


;plot, restframespectrum[0,*], restframespectrum[1,*]

specphot, restframespectrum, restframemags



; check to see if the spectrum is a string (ie a filename) or an array
s=size(spectrum)
if (s[1] eq 7) then begin
	readcol,restframespectrum,sp_wave,sp_flux,/silent
endif else begin
	sp_wave=restframespectrum[0,*]
	sp_flux=restframespectrum[1,*]
endelse

sp_wave_z=sp_wave/(1+z)

;obsframespectrum=[ transpose(sp_wave_z), transpose(sp_flux)]
obsframespectrum=[ (sp_wave_z), (sp_flux)]

;oplot, obsframespectrum[0,*], obsframespectrum[1,*]


specphot, obsframespectrum, obsframemags

	kcor=obsframemags-restframemags+2.5*alog10(1+z)
end

;;;;;;;;;;;;;;;;;
pro fluxfactors, spectrum, fluxfactors

h = 1D*6.626E-27
c = 1D*2.998E18
hc = 1D*1.986E-8
e = 1D*2.71828
k = 1D*1.38E-16

;;Read in the filter Effective Area curves
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/V_UVOT.txt",lambda,V_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/B_UVOT.txt",lambda,B_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/U_UVOT.txt",lambda,U_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVW1.txt",lambda,W1_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVM2.txt",lambda,M2_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVW2.txt",lambda,W2_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVW1-rc.txt",lambda,W1rc_EA,/silent
readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/UVW2-rc.txt",lambda,W2rc_EA,/silent

; check to see if it is a string (ie a filename) or an array
s=size(spectrum)
if (s[1] eq 7) then begin
	readcol,spectrum,sp_wave,sp_flux,/silent
endif else begin
	sp_wave=spectrum[0,*]
	sp_flux=spectrum[1,*]
endelse

flux=interpol(sp_flux,sp_wave,lambda)

v_counts=V_EA*flux*(10*lambda/hc)
b_counts=B_EA*flux*(10*lambda/hc)
u_counts=U_EA*flux*(10*lambda/hc)
w1_counts=W1_EA*flux*(10*lambda/hc)
m2_counts=M2_EA*flux*(10*lambda/hc)
w2_counts=W2_EA*flux*(10*lambda/hc)
w1rc_counts=W1rc_EA*flux*(10*lambda/hc)
w2rc_counts=W2rc_EA*flux*(10*lambda/hc)

v_spec=-2.5*alog10(total(v_counts))+17.89
b_spec=-2.5*alog10(total(b_counts))+19.11
u_spec=-2.5*alog10(total(u_counts))+18.34
w1_spec=-2.5*alog10(total(w1_counts))+17.49
m2_spec=-2.5*alog10(total(m2_counts))+16.82
w2_spec=-2.5*alog10(total(w2_counts))+17.35
w1rc_spec=-2.5*alog10(total(w1rc_counts))+17.34
w2rc_spec=-2.5*alog10(total(w2rc_counts))+17.25

vband=where(lambda ge 5400 and lambda le 5500)
bband=where(lambda ge 4300 and lambda le 4400)
uband=where(lambda ge 3400 and lambda le 3500)
w1band=where(lambda ge 2550 and lambda le 2650)
m2band=where(lambda ge 2200 and lambda le 2300)
w2band=where(lambda ge 1900 and lambda le 2000)

counts=[total(w2_counts),total(m2_counts),total(w1_counts),total(u_counts),total(b_counts),total(v_counts)]

flatflux=[mean(flux[w2band]),mean(flux[m2band]),mean(flux[w1band]),mean(flux[uband]),mean(flux[bband]),mean(flux[vband])]

fluxfactors=flatflux/counts

filterflux=[2.728, 1.391, 1.666, 4.027, 9.179, 7.024]*10.0^(-16.0)
filterflux=[7.024, 9.179, 4.027, 1.666, 1.391, 2.728]*10.0^(-16.0)


;print, fluxfactors

;print, "compare to nominal"

;print, filterflux


;sedflat=[ [1900,flux[5]], [2000, flux[5] ], [2200, flux[4] ], [2300, flux[4] ], [2550,flux[3]], [2650,flux[3]], $ 
;	[3400, flux[2]], [3500, flux[2]], [4300, flux[1]], [4400,flux[1]], [5400,flux[0]],  [5500,flux[0] ] ]

;plot, sedflat[0,*], sedflat[1,*]

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hsiao_all

z=0
;distmod=29.75
;distmoderr=0.12
ebvmw=0.00
ebvmwerr=0.00
ebvhost=0.0
ebvhosterr=0.00


readcol,"/uufs/astro.utah.edu/common/home/u0676985/SN/hsiao/snflux_1a.dat",epochs, sp_wave, sp_flux,/silent

JD=fltarr(105)

mag_array=fltarr(8,105)
kcor_array=fltarr(8,n_elements(JD))
redmag_array=fltarr(2,n_elements(JD))
redfrac_array=fltarr(2,n_elements(JD))
snspecmags_array=fltarr(8,n_elements(JD))
RlambdaMW=fltarr(8,n_elements(JD))
RlambdaCCLMC=fltarr(8,n_elements(JD))
fluxfactors_array=fltarr(6,n_elements(JD))



for i=0,104 do begin
epoch=i-20
JD(i)=epoch
range=where(epochs eq epoch)


spec=[transpose(sp_wave[range]),transpose(sp_flux[range])]

	specphot, spec, mags
	mag_array[*,i]=mags

;readcol,uvotfile,JD, w2, m2, w1, uu, bb, vv, /silent

readcol,"/uufs/astro.utah.edu/common/home/u0676985/snscripts/V_UVOT.txt",lambda,V_EA,/silent

h = 1D*6.626E-27
c = 1D*2.998E18
hc = 1D*1.986E-8
e = 1D*2.71828
k = 1D*1.38E-16

snmags=mag_array(0:5,*)

;; compute input spectrum based on photometry

;sed, snmags[*,i], sedflat=sedflat

;spec=sedflat


spec_file=spec

flux=interpol(sp_flux,sp_wave,lambda)

Ahost=sne_goobarlmc_reddening(lambda,Ebvhost)
spechostdred=flux*10^(-Ahost/2.5)
; for flat spechostdred=sp_flux*0+1*10^(-Ahost/2.5)

;  redshift
lambda_red=lambda*(1+z)

Amw=sne_mw_reddening(lambda_red,Ebvmw)
;Amw=sne_mw_reddening(sp_wave_red,0)
specmwdred=spechostdred*10^(-Amw/2.5)

spec=fltarr(2,n_elements(lambda))
for j=0,n_elements(lambda)-1 do spec[0,j]=lambda[j]
for j=0,n_elements(lambda)-1 do spec[1,j]=specmwdred[j]

plot, lambda, specmwdred

mangle, spec, snmags[*,i], z, mangledspectrum

kcor, mangledspectrum, z, kcor
kcor_array[*,i]=kcor

specphot, mangledspectrum, mangledmags, redfrac=redfrac

redmag_array[*,i]=[mangledmags[2]-mangledmags[6], mangledmags[0]-mangledmags[7]]
redfrac_array[*,i]=redfrac

fluxfactors, mangledspectrum, fluxfactors
fluxfactors_array[*,i]=fluxfactors

snspecmags_array[*,i]=mangledmags

;;;; compute extinction coefficients here
Amw=sne_mw_reddening(mangledspectrum[0,*],0.1)
specmwdred=mangledspectrum[1,*]*10^(-Amw/2.5)

specphot, [ mangledspectrum[0,*],(specmwdred)], deredmags

RlambdaMW[*,i]= (deredmags - mangledmags)/0.1

Ahost=sne_goobarlmc_reddening(mangledspectrum[0,*],0.1)
spechostdred=mangledspectrum[1,*]*10^(-Ahost/2.5)

specphot, [ mangledspectrum[0,*],(spechostdred)], deredmags

RlambdaCCLMC[*,i]= (deredmags - mangledmags)/0.1



;stop

endfor


SET_PLOT, 'PS'
DEVICE, FILENAME='hsiao_all.ps',/LANDSCAPE, /COLOR
;!P.MULTI = [4,2,2]
orchid  = FSC_COLOR('Orchid',!D.TABLE_SIZE-2)
blue  = FSC_COLOR('Blue',!D.TABLE_SIZE-3)
cyan  = FSC_COLOR('Cyan',!D.TABLE_SIZE-4)
green  = FSC_COLOR('Green',!D.TABLE_SIZE-5)
orange  = FSC_COLOR('Orange',!D.TABLE_SIZE-6)
red  = FSC_COLOR('Red',!D.TABLE_SIZE-7)

allsize=1.25
vvsym=-4
bbsym=-5
uusym=-6
w1sym=-7
m2sym=-4
w2sym=-5
w1rcsym=w1sym
w2rcsym=w2sym
w1rcsize=allsize+0.5
w2rcsize=allsize+0.5
vvcolor=green
bbcolor=blue
uucolor=orchid
w1color=cyan
m2color=red
w2color=orange
w1rccolor=0
w2rccolor=0

plot, JD, snmags[5,*], PSYM=vvsym, THICK=1,$
	XTHICK=3, YTHICK=3, CHARTHICK=3, SYMSIZE=allsize, /XSTYLE, /YSTYLE,$
	XRANGE=[-20,85], YRANGE=[2,-8], XTITLE='Epoch',$
	YTITLE='Vega Magnitudes', CHARSIZE=1.5;, COLOR=vvcolor
oplot, JD, snmags[4,*], PSYM=bbsym, SYMSIZE=allsize, COLOR=bbcolor
oplot, JD, snmags[3,*], PSYM=uusym, SYMSIZE=allsize, COLOR=uucolor
oplot, JD, snmags[2,*], PSYM=w1sym, SYMSIZE=allsize, COLOR=w1color
oplot, JD, snmags[1,*], PSYM=m2sym, SYMSIZE=allsize;, COLOR=m2color
oplot, JD, snmags[0,*], PSYM=w2sym, SYMSIZE=allsize;, COLOR=w2color
oplot, JD, mag_array[6,*], PSYM=w1rcsym, SYMSIZE=allsize;, COLOR=m2color
oplot, JD, mag_array[7,*], PSYM=w2rcsym, SYMSIZE=allsize;, COLOR=w2color


plot, JD, kcor_array[5,*], PSYM=vvsym, THICK=1,$
	XTHICK=3, YTHICK=3, CHARTHICK=3, SYMSIZE=allsize, /XSTYLE, /YSTYLE,$
	XRANGE=[-20,85], YRANGE=[-0.02,0.02], XTITLE='Epoch',$
	YTITLE='k corrections', CHARSIZE=1.5;, COLOR=vvcolor
oplot, JD, kcor_array[4,*], PSYM=bbsym, SYMSIZE=allsize, COLOR=bbcolor
oplot, JD, kcor_array[3,*], PSYM=uusym, SYMSIZE=allsize, COLOR=uucolor
oplot, JD, kcor_array[2,*], PSYM=w1sym, SYMSIZE=allsize, COLOR=w1color
oplot, JD, kcor_array[1,*], PSYM=m2sym, SYMSIZE=allsize, COLOR=m2color
oplot, JD, kcor_array[0,*], PSYM=w2sym, SYMSIZE=allsize, COLOR=w2color
oplot, JD, kcor_array[6,*], PSYM=w1rcsym, SYMSIZE=w1rcsize;, COLOR=w1rccolor
oplot, JD, kcor_array[7,*], PSYM=w2rcsym, SYMSIZE=w2rcsize;, COLOR=w2rccolor



plot, JD, RlambdaMW[5,*], PSYM=vvsym, THICK=1,$
	XTHICK=3, YTHICK=3, CHARTHICK=3, SYMSIZE=allsize, /XSTYLE, /YSTYLE,$
	XRANGE=[-20,85], YRANGE=[0,16], XTITLE='Epoch',$
	YTITLE='Rlambda', CHARSIZE=1.5;, COLOR=vvcolor
oplot, JD, RlambdaMW[4,*], PSYM=bbsym, SYMSIZE=allsize, COLOR=bbcolor
oplot, JD, RlambdaMW[3,*], PSYM=uusym, SYMSIZE=allsize, COLOR=uucolor
oplot, JD, RlambdaMW[2,*], PSYM=w1sym, SYMSIZE=allsize, COLOR=w1color
oplot, JD, RlambdaMW[1,*], PSYM=m2sym, SYMSIZE=allsize, COLOR=m2color
oplot, JD, RlambdaMW[0,*], PSYM=w2sym, SYMSIZE=allsize, COLOR=w2color
oplot, JD, RlambdaMW[6,*], PSYM=w1rcsym, SYMSIZE=w1rcsize;, COLOR=w1rccolor
oplot, JD, RlambdaMW[7,*], PSYM=w2rcsym, SYMSIZE=w2rcsize;, COLOR=w2rccolor


plot, JD, RlambdaCCLMC[5,*], PSYM=vvsym, THICK=1,$
	XTHICK=3, YTHICK=3, CHARTHICK=3, SYMSIZE=allsize, /XSTYLE, /YSTYLE,$
	XRANGE=[-20,85], YRANGE=[0,16], XTITLE='Epoch',$
	YTITLE='Rlambda', CHARSIZE=1.5;, COLOR=vvcolor
oplot, JD, RlambdaCCLMC[4,*], PSYM=bbsym, SYMSIZE=allsize, COLOR=bbcolor
oplot, JD, RlambdaCCLMC[3,*], PSYM=uusym, SYMSIZE=allsize, COLOR=uucolor
oplot, JD, RlambdaCCLMC[2,*], PSYM=w1sym, SYMSIZE=allsize, COLOR=w1color
oplot, JD, RlambdaCCLMC[1,*], PSYM=m2sym, SYMSIZE=allsize, COLOR=m2color
oplot, JD, RlambdaCCLMC[0,*], PSYM=w2sym, SYMSIZE=allsize, COLOR=w2color
oplot, JD, RlambdaCCLMC[6,*], PSYM=w1rcsym, SYMSIZE=w1rcsize;, COLOR=w1rccolor
oplot, JD, RlambdaCCLMC[7,*], PSYM=w2rcsym, SYMSIZE=w2rcsize;, COLOR=w2rccolor

fluxfactors_array=fluxfactors_array*10^16

plot, JD, fluxfactors_array[5,*], PSYM=vvsym, THICK=1,$
	XTHICK=3, YTHICK=3, CHARTHICK=3, SYMSIZE=allsize, /XSTYLE, /YSTYLE,$
	XRANGE=[-20,85], YRANGE=[0,1], XTITLE='Epoch',$
	YTITLE='flux conversion factors', CHARSIZE=1.5;, COLOR=vvcolor
oplot, JD, fluxfactors_array[4,*], PSYM=bbsym, SYMSIZE=allsize, COLOR=bbcolor
oplot, JD, fluxfactors_array[3,*], PSYM=uusym, SYMSIZE=allsize, COLOR=uucolor
oplot, JD, fluxfactors_array[2,*], PSYM=w1sym, SYMSIZE=allsize, COLOR=w1color
oplot, JD, fluxfactors_array[1,*], PSYM=m2sym, SYMSIZE=allsize, COLOR=m2color
oplot, JD, fluxfactors_array[0,*], PSYM=w2sym, SYMSIZE=allsize, COLOR=w2color

plot, JD, redfrac_array[0,*], PSYM=w1rcsym, THICK=1,$
	XTHICK=3, YTHICK=3, CHARTHICK=3, SYMSIZE=w1rcsize, /XSTYLE, /YSTYLE,$
	XRANGE=[-20,85], YRANGE=[0,1], XTITLE='Epoch',$
	YTITLE='UV fraction', CHARSIZE=1.5, COLOR=w1rccolor
oplot, JD, redfrac_array[1,*], PSYM=w2rcsym, SYMSIZE=w2rcsize, COLOR=w2rccolor

leg=['vv', 'bbs', 'uu', 'w1', 'm2', 'w2', 'w1rc', 'w2rc']
legend, leg, psym=[vvsym, bbsym, uusym, w1sym, m2sym, w2sym, w1rcsym, w2rcsym], $
symsize=[allsize, allsize, allsize, allsize, allsize, allsize, w1rcsize, w2rcsize], $
	color=[vvcolor, bbcolor, uucolor, w1color, m2color, w2color, w1rccolor, w2rccolor], pos=[0.5,0.45], /norm, charsize=1.5


DEVICE, /CLOSE
SET_PLOT, 'X'
 $gv hsiao_all.ps &

save, /all, filename='hsiao_all.sav'


print, "last stop"
;stop
end

