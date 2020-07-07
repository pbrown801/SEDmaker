function zpphot, spectrum, filter, zp

;;Read in the filter Effective Area curves
;readcol,"$SNSCRIPTS/V_UVOT.txt",lambda,V_EA,/silent
;readcol,"$SNSCRIPTS/B_UVOT.txt",lambda,B_EA,/silent
;readcol,"$SNSCRIPTS/U_UVOT.txt",lambda,U_EA,/silent
;readcol,"$SNSCRIPTS/UVW1.txt",lambda,W1_EA,/silent
;readcol,"$SNSCRIPTS/UVM2.txt",lambda,M2_EA,/silent
;readcol,"$SNSCRIPTS/UVW2.txt",lambda,W2_EA,/silent
;readcol,"$SNSCRIPTS/UVW1-rc.txt",lambda,W1rc_EA,/silent
;readcol,"$SNSCRIPTS/UVW2-rc.txt",lambda,W2rc_EA,/silent
; 
; v_spec=-2.5*alog10(total(v_counts))+17.89
; b_spec=-2.5*alog10(total(b_counts))+19.11
; u_spec=-2.5*alog10(total(u_counts))+18.34
; w1_spec=-2.5*alog10(total(w1_counts))+17.49
; m2_spec=-2.5*alog10(total(m2_counts))+16.82
; w2_spec=-2.5*alog10(total(w2_counts))+17.35
; w1rc_spec=-2.5*alog10(total(w1rc_counts))+17.34
; w2rc_spec=-2.5*alog10(total(w2rc_counts))+17.25



; $SNSCRIPTS/csp/CSP_filter_package/V_LC9844_texas_WLcorr.txt
;
;
;u 3639.3 13.044
;B 4350.6 14.344
;V-3014 5369.6 14.540
;V-3099 5401.4 	14.493
;V-9844 5375.2 14.450
;g 4765.1 15.135
;r 6223.3 14.915
;i 7609.2 14.781
;Y 10289.1 12.687
;J 12382.0 12.853
;H 16272.2 12.555
;K 21422.2 11.967

h = 1D*6.626E-27
c = 1D*2.998E18
hc = 1D*1.986E-8
e = 1D*2.71828
k = 1D*1.38E-16

;;;;;;;;;; find the vega zeropoint for the filter
;;Read in the filter  curves
readcol, filter, filterlambda,filtertransmission,/silent

;;;;; put the filter transmission curve on a 10 A grid

wavestart=floor(min(filterlambda,/nan)/10)*10
waveend=floor(max(filterlambda,/nan)/10+1)*10

lambda=intarr( (waveend-wavestart)/10)
for i=0,n_elements(lambda)-1 do lambda[i]=wavestart+10*i

transmission=interpol(filtertransmission,filterlambda,lambda)

;; Set the zero point using vega
;;Read in the Spectra and interpolate to the 10 A grid

;readcol,'$SNSCRIPTS/vega.dat',sp_wave,sp_flux,/silent
;zero=where(finite(sp_flux) eq 0)
;sp_flux[zero]=0.0

;vegaflux=interpol(sp_flux,sp_wave,lambda)
;vega_counts=transmission*vegaflux*(10.0*lambda/hc)
;vega_mag=-2.5*alog10(total(vega_counts))


;;;;;;;;; add errors if spectrum doesn't cover filter or vice versa

;;;;;;;;;;;;;;;;;;

;;; compute magnitude from source spectrum in the above filter

; check to see if it is a string (ie a filename) or an array
s=size(spectrum)
if (s[1] eq 7) then begin
	readcol,spectrum,sp_wave,sp_flux,/silent
endif else begin
	sp_wave=spectrum[0,*]
	sp_flux=spectrum[1,*]
endelse
zero=where(finite(sp_flux) eq 0)
sp_flux[zero]=0.0

flux=interpol(sp_flux,sp_wave,lambda)
filter_counts=transmission*flux*(10.0*lambda/hc)

;filter_mag=-2.5*alog10(total(filter_counts))-vega_mag
filter_mag=-2.5*alog10(total(filter_counts))+zp
return, filter_mag


end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hsiao_scor

readcol,"snflux_1a.dat",epochs, sp_wave, sp_flux,/silent

epoch_array=fltarr(105)

v_uvot_array=fltarr(105)
v_uvot_redo_array=fltarr(105)
v_uvot_flatter_array=fltarr(105)
v_csp_array=fltarr(105)
v_csp_noatm_array=fltarr(105)

for i=0,104 do begin
	epoch=i-20
	epoch_array(i)=epoch
	range=where(epochs eq epoch)
	spec=[transpose(sp_wave[range]),transpose(sp_flux[range])]

	v_uvot=zpphot(spec,"$SNSCRIPTS/V_UVOT.txt", 17.89)
	v_uvot_array[i]=v_uvot

	v_uvot_redo=zpphot(spec,"$SNSCRIPTS/V_UVOT_redo.txt", 18.2775)
	v_uvot_redo_array[i]=v_uvot_redo

	v_uvot_flatter=zpphot(spec,"$SNSCRIPTS/V_UVOT_flatter.txt", 19.5891)
	v_uvot_flatter_array[i]=v_uvot_flatter

	v_csp=zpphot(spec,'$SNSCRIPTS/csp/CSP_filter_package/V_LC9844_texax_WLcorr_atm.txt', 14.450)
	v_csp_array[i]=v_csp

	v_csp_noatm=zpphot(spec,'$SNSCRIPTS/csp/CSP_filter_package/V_LC9844_texas_WLcorr.txt', 14.637)
	v_csp_noatm_array[i]=v_csp_noatm

endfor

readcol,"$SNSCRIPTS/V_UVOT.txt",lambda,V_EA,/silent
readcol,"$SNSCRIPTS/V_UVOT_redo.txt",lambda,V_redo_EA,/silent


stop
SET_PLOT, 'PS'
DEVICE, FILENAME='hsiao_scor.ps',/LANDSCAPE, /COLOR

plot, epoch_array, v_csp_array-v_uvot_array, PSYM=4, THICK=1,$
	XTHICK=3, YTHICK=3, CHARTHICK=3, SYMSIZE=allsize, /XSTYLE, /YSTYLE,$
	XRANGE=[-20,85], YRANGE=[-0.05,0.15], XTITLE='Epoch',$
	YTITLE='v_csp - v_uvot', CHARSIZE=1.5;, COLOR=vvcolor

plot, epoch_array, v_uvot_array-v_uvot_redo_array, PSYM=4, THICK=1,$
	XTHICK=3, YTHICK=3, CHARTHICK=3, SYMSIZE=allsize, /XSTYLE, /YSTYLE,$
	XRANGE=[-20,85], YRANGE=[-0.05,0.15], XTITLE='Epoch',$
	YTITLE='v_uvot - v_uvot_redo', CHARSIZE=1.5;, COLOR=vvcolor

plot, epoch_array, v_uvot_array-v_uvot_flatter_array, PSYM=4, THICK=1,$
	XTHICK=3, YTHICK=3, CHARTHICK=3, SYMSIZE=allsize, /XSTYLE, /YSTYLE,$
	XRANGE=[-20,85], YRANGE=[-0.05,0.15], XTITLE='Epoch',$
	YTITLE='v_uvot - v_uvot_flatter', CHARSIZE=1.5;, COLOR=vvcolor

plot, epoch_array, v_csp_array-v_csp_noatm_array, PSYM=4, THICK=1,$
	XTHICK=3, YTHICK=3, CHARTHICK=3, SYMSIZE=allsize, /XSTYLE, /YSTYLE,$
	XRANGE=[-20,85], YRANGE=[-0.05,0.15], XTITLE='Epoch',$
	YTITLE='v_csp - v_csp_noatm', CHARSIZE=1.5;, COLOR=vvcolor



DEVICE, /CLOSE
SET_PLOT, 'X'
 $open hsiao_scor.ps &

save, /all, filename='hsiao_scor.sav'


print, "last stop"
stop
end

