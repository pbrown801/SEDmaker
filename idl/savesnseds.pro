pro savesnseds



SNlist=['ASASSN-15lh','SN2011by','SN2007af', 'SN2005ke', 'SN2005hk', 'SN2011de','SN2007Y',  'SN2006jc',   'SN2006aj',  'SN2011dh', 'SN2012aw', 'SN2008aw',   'SN2010jl', 'SN2008es']
nSNe=n_elements(SNlist)

SNlegend=['ASASSN-15lh - SLSN I', 'SN2011by - Ia-blue','SN2007af - Ia-red', 'SN2005ke - Ia-91bg', 'SN2012Z - Ia-02cx','SN2011de - Ia SC','SN2007Y - Ib',  'SN2006jc - Ibn',  'SN2006aj - Ic/GRB',  'SN2011dh - IIb', 'SN2012aw - IIP', 'SN2008aw - IIL',   'SN2010jl - IIn', 'SN2008es - SLSNII']
SNexplosiondates=[55796.687,54157.12, 53685.77, 53685.1-15.0,55695,54163.3-20.0,  54018,  53784.1,  55712.5, 56002.0, 54526.5,   55470+30, 54590]
SNexplosionrefs=['Nugent_etal_2011','Brown_etal_2012a','Brown_etal_2012a', 'Phillips_etal_2007','Brown et al. 2013 peak -20','Stritzinger_etal_2009 B peak - 20', 'Nakana_etal_2006 a few days before discovery',   'Campana_etal_2006 GRB trigger time',  'Arcavi_etal_2011', 'Fraser_etal_2012', 'discover - 1.5 days',   'thirty days before first observation Stoll et al. 2011', 'ten days before first observation']

colors=['red','maroon', 'orange', 'orange red','magenta','dark green', 'lawn green',  'black', 'turquoise', 'blue', 'royal blue',   'purple', 'grey']

symnum=[14,15,16,4,5,6,7,9,11]
sizes=[1,1,1,1,1,1,1,1,1] ;,1,1,1]

nprimeSNe=3

colors=['blue','red','green','red','red','blue','blue','green','blue','red','green', $
'blue','red','green','blue','red','green','blue','red','green','blue','red','green' ] 
primecolors=['black', 'black','black' ] 
primecolors=['blue', 'red','dark green' ] 
primecolors=primecolors[0:nprimeSNe-1]
colors=colors[nprimeSNe:nSNe+nprimeSNe-1]
allcolors=[primecolors,colors]

altcolors=['powder blue', 'pink', 'light sea green', 'powder blue', 'pink', 'light sea green', 'powder blue', 'pink', 'light sea green', 'powder blue', 'pink', 'light sea green', 'powder blue', 'pink', 'light sea green', 'powder blue', 'pink', 'light sea green', 'powder blue', 'pink', 'light sea green', 'powder blue', 'pink', 'light sea green']

primesymnum=[14,15,16,14,15,16]
symnum=[4,5,6,9,7,11,23,25,27,29,35,37,41,4,5,6,9,7,11,23,25,27,29,35,37,41]
symnum=symnum[0:nSNe-1]
primesymnum=primesymnum[0:nprimeSNe-1]
allsymnum=[primesymnum,symnum]

primesizes=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
sizes=[0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8]
primesizes=primesizes[0:nprimeSNe-1]
sizes=sizes[0:nSNe-1]
allsizes=[primesizes,sizes]


sizes=sizes[0:nSNe-1]




;;;;;

SNlist=['ASASSN-15lh','SN2011by','SN2007af', 'SN2005ke', 'SN2005hk', 'SN2011de','SN2007Y',  'SN2006jc',   'SN2006aj',  'SN2011dh', 'SN2012aw', 'SN2008aw',   'SN2010jl', 'SN2008es']


spectralist=['ASASSN-15lh','SN2011by','SN2007af', 'SN2005ke', 'SN2005hk', 'SN2011aa','SN2007Y',  'SN2006jc',   'SN2006aj',  'SN2011dh', 'SN2012aw', 'SN2008aw',   'SN2010jl', 'SN2008es']



spectralist2=['$SNSCRIPTS/vega.dat','$SNSCRIPTS/SN1992A_UV.dat','$SNSCRIPTS/SN1999em_UV.dat']


SNlegend=['ASASSN-15lh - SLSN I','SN2011by - Ia-blue','SN2007af - Ia-red', 'SN2005ke - Ia-91bg', 'SN2012Z - Ia-02cx','SN2011de - Ia SC','SN2007Y - Ib',  'SN2006jc - Ibn',  'SN2006aj - Ic/GRB',  'SN2011dh - IIb', 'SN2012aw - IIP', 'SN2008aw - IIL',   'SN2010jl - IIn', 'SN2008es - SLSN II']
SNexplosiondates=[57190.0,55675,54157.12, 53685.77, 53685.1-15.0,55695,54163.3-20.0,  54018,  53784.1,  55712.5, 56002.0, 54526.5,   55470+30, 54590]
SNexplosionrefs=['picked to fit', 'picked to fit','Brown_etal_2012a','Brown_etal_2012a', 'Phillips_etal_2007','Brown et al. 2013 peak -20','Stritzinger_etal_2009 B peak - 20', 'Nakana_etal_2006 a few days before discovery',   'Campana_etal_2006 GRB trigger time',  'Arcavi_etal_2011', 'Fraser_etal_2012', 'discover - 1.5 days',   'thirty days before first observation Stoll et al. 2011', 'ten days before first observation']

colors=['black','red','maroon', 'orange', 'orange red','magenta','dark green', 'lawn green',  'black', 'turquoise', 'blue', 'royal blue',   'purple', 'grey']

symnum=[14,15,16,4,5,6,7,9,11]
sizes=[1,1,1,1,1,1,1,1,1] ;,1,1,1]

nSNe=n_elements(SNlist)
nprimeSNe=3

colors=['blue','red','green','red','red','blue','blue','green','blue','red','green', $
'blue','red','green','blue','red','green','blue','red','green','blue','red','green' ] 
primecolors=['black', 'black','black' ] 
primecolors=['black','blue', 'red','dark green' ] 
primecolors=primecolors[0:nprimeSNe-1]
colors=colors[nprimeSNe:nSNe+nprimeSNe-1]
allcolors=[primecolors,colors]

altcolors=['black','powder blue', 'pink', 'light sea green', 'powder blue', 'pink', 'light sea green', 'powder blue', 'pink', 'light sea green', 'powder blue', 'pink', 'light sea green', 'powder blue', 'pink', 'light sea green', 'powder blue', 'pink', 'light sea green', 'powder blue', 'pink', 'light sea green', 'powder blue', 'pink', 'light sea green']

primesymnum=[14,15,16,14,15,16]
symnum=[4,5,6,9,7,11,23,25,27,29,35,37,41,4,5,6,9,7,11,23,25,27,29,35,37,41]
symnum=symnum[0:nSNe-1]
primesymnum=primesymnum[0:nprimeSNe-1]
allsymnum=[primesymnum,symnum]

primesizes=[1.5,1,1,1,1,1,1,1,1,1,1,1,1,1,11,1,1]
sizes=[0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8]
primesizes=primesizes[0:nprimeSNe-1]
sizes=sizes[0:nSNe-1]
allsizes=[primesizes,sizes]


sizes=sizes[0:nSNe-1]



nSNe=n_elements(SNlist)

nepochs_array=intarr(nSNe)

for s=0, nSNe-1 do begin
SNname=SNlist[s]
restore, filename=SNname+'_sed.sav'

nallepochs=n_elements( dt.time_array[*] )
goodfilters=make_array(nallepochs, /float, value=!Values.F_NAN)

for n=0,nallepochs-1 do goodfilters[n]=total(finite(dt.counts_array[*,n]))
goodepochs=where(goodfilters eq 6)
nepochs=n_elements(goodepochs)
sedmjd_array=dt.time_array[goodepochs]

nepochs_array[s]=nepochs

endfor

mostepochs=max(nepochs_array)


allsnflux_array=make_array(nSNe,4,mostepochs, /float, value=!Values.F_NAN)
allsnluminosity_array=make_array(nSNe,mostepochs, /float, value=!Values.F_NAN)
allsnlumlog_array=make_array(nSNe,mostepochs, /float, value=!Values.F_NAN)


for s=0, nSNe-1 do begin
	SNname=SNlist[s]
	restore, filename=SNname+'_sed.sav'

	nallepochs=n_elements( dt.time_array[*] )
	goodfilters=make_array(nallepochs, /float, value=!Values.F_NAN)

	for n=0,nallepochs-1 do goodfilters[n]=total(finite(dt.counts_array[*,n]))
	goodepochs=where(goodfilters eq 6)
	nepochs=n_elements(goodepochs)
	sedmjd_array=dt.time_array[goodepochs]

	nepochs_array[s]=nepochs

	allsnflux_array[s,0,0:nepochs-1]=sedmjd_array
	allsnflux_array[s,1,0:nepochs-1]=sedarray.muvflux_array
	allsnflux_array[s,2,0:nepochs-1]=sedarray.nuvflux_array
	allsnflux_array[s,3,0:nepochs-1]=sedarray.optflux_array
	allsnluminosity_array[s,0:nepochs-1]=sedarray.luminosity
	allsnlumlog_array[s,0:nepochs-1]=sedarray.luminosity_log

endfor

save, 'snseds.sav', allsnflux_array, allsnluminosity_array, allsnlumlog_array, snexplosiondates



print, 'final stop'
stop
end
