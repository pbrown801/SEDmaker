pro sbck_makespec

allfiles='sbck_spectratrimmed.txt'

readcol,allfiles, fitsspectra, format='A'

for j=0,n_elements(fitsspectra) - 1 do makespecfromfits, fitsspectra[j]

stop
end