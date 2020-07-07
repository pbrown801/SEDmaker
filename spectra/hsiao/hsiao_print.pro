;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hsiao_print

z=0
;distmod=29.75
;distmoderr=0.12
ebvmw=0.00
ebvmwerr=0.00
ebvhost=0.0
ebvhosterr=0.00


readcol,"$SNFOLDER/hsiao/snflux_1a.dat",epochs, sp_wave, sp_flux,/silent

JD=intarr(105)

for i=0,104 do begin
	epoch=i-20
	JD(i)=epoch
	range=where(epochs eq epoch)

	filename='hsiaospec/hsiao_'+strtrim(string(epoch),1)+'.dat'


	openw, lun2, filename, /get_lun

	for q=0,n_elements(range)-1 do 	printf, lun2, sp_wave[range[q]],sp_flux[range[q]]

	close, lun2
	free_lun, lun2




endfor




stop
end

