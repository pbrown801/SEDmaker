;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro scaletemplate

trange=120.0

filename='snflux_1a.dat'

readcol, filename,epochs, sp_wave, sp_flux,/silent

n_epochs=n_elements(where(epochs[uniq(epochs)] lt trange))

epoch_array=epochs[uniq(epochs)]
epoch_array=epoch_array[0:n_epochs-1]


maxflux=max(sp_flux)

openw,lun1, 'snflux_1a_scaled.txt', /get_lun

for i=0,n_epochs-1 do begin

	range=where(epochs eq epoch_array[i] and sp_wave le 11000.0 and sp_wave ge 1000.0 )
	; some epochs are duplicated
	for r=0, n_elements(range)-1 do printf, lun1, (epochs[range[r]]+20.0)*(10.0/trange), ' , ', (sp_wave[range[r]]-1000)/10000*10.0, ' , ', sp_flux[range[r]]*100.0/maxflux
print, epoch_array[i]

endfor

close, lun1
free_lun, lun1


stop
end

