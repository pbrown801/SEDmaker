pro makesnseds

SNlist=['SN2011by','SN2007af', 'SN2005ke', 'SN2005hk', 'SN2011de','SN2007Y',  'SN2006jc',   'SN2006aj',  'SN2011dh', 'SN2012aw', 'SN2008aw',   'SN2010jl', 'SN2008es']

SNlist=['SN2005cf','SN2006ej','SN2007af','SN2007co', 'SN2008Q', 'SN2008ec',  'SN2008hv',  'SNF20080514002']

SNlist=['SN2012hr',  'SN2011ia',  'SN2011fe',  'SN2011dn',  'SN2011at',  'SN2011ao',  'SN2011B',  'SN2010ev', 'SN2009cz', 'SN2006dm']


for n=0, n_elements(SNlist) -1 do makesedarray15, SNlist[n], sedarray



print, 'final stop'
stop
end
