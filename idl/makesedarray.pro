pro makesedarray, SNname, sedarray


photfile='$SOUSA/data/'+SNname+'_uvotB14.1.dat'
pjb_phot_array_B141, photfile,   dt=dt

nallepochs=n_elements( dt.time_array[*] )
goodfilters=make_array(nallepochs, /float, value=!Values.F_NAN)

for n=0,nallepochs-1 do goodfilters[n]=total(finite(dt.counts_array[*,n]))
goodepochs=where(goodfilters eq 6)
nepochs=n_elements(goodepochs)
muvflux_array=make_array(nepochs, /float, value=!Values.F_NAN)
nuvflux_array=make_array(nepochs, /float, value=!Values.F_NAN)
optflux_array=make_array(nepochs, /float, value=!Values.F_NAN)
sed_array=make_array(nepochs, 7,/float, value=!Values.F_NAN)
sedmjd_array=dt.time_array[goodepochs]
for n=0,nepochs-1 do begin

	inputcounts=dt.counts_array[*,goodepochs[n]]
	inputcounterrs=dt.countserr_array[*,goodepochs[n]]

		gridsed, inputcounts, inputcounterrs=inputcounterrs, bestsed=bestsed, sedarray=sedarray, chisqarray=chisqarray, sedcount_array=sedcount_array, sigmaone=sigmaone, sedwave=sedwave

		sed_array[n,*]=bestsed

		uvfluxfraction, [transpose(sedwave),transpose(bestsed)], muvflux, nuvflux, optflux

		muvflux_array[n]=muvflux
		nuvflux_array[n]=nuvflux
		optflux_array[n]=optflux

endfor

restore, '$SNFOLDER/www/SwiftSN/host.sav'

index=where(host.snname_array eq SNname)
mu_best=host.dm_best_array[index]

distance_parsec=10.0^((mu_best/5.0)+1.0)
;;; distance modulus to centimeters/10^18
distance_cm=1.0d0
distance_cm=distance_parsec*3.08568*10.0^18.0
distance_cm18=distance_parsec*3.08568
distance=distance_parsec*3.08568

luminosity=distance^2.0*4.0*3.14159265/10.0^6.0*(muvflux_array+optflux_array+nuvflux_array)
luminosity_36=distance_cm18[0]^2.0*4.0*3.14159265*(muvflux_array+optflux_array+nuvflux_array)
luminosity_log=alog10(luminosity_36)+36.0

sedarray = CREATE_STRUCT('SNname', SNname, $
'mjd_array', dt.time_array, $
'mag_array', dt.mag_array, $
'counts_array', dt.counts_array, $
'countserr_array', dt.countserr_array, $
'sedmjd_array', sedmjd_array, $
'sed_array', sed_array, $
'sedwave', sedwave, $
'muvflux_array', muvflux_array, $
'nuvflux_array', nuvflux_array, $
'optflux_array', optflux_array, $
'distance_parsec', distance_parsec, $
'mu_best', mu_best, $
'distance_cm', distance_cm, $
'luminosity_log', luminosity_log, $
'luminosity', luminosity )


save, filename=SNname+'_sed.sav', SNname, sedarray, dt

end
