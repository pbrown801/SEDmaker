pro plotsnseds



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
xrange=[-5,50]
yrangetop=[0.0,0.55]
yrangebottom=[0.0,0.12]

xdata=[1,2,3,4]
ydata=[2,3,4,5]

nxticks=11
nyticks=3

figurename='boloflux3plots_revised.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=12, bits_per_pixel=8, /color


q=0


;;;;;;;;;;;;;;;;;;;;;;;;;

plot, xdata, ydata, /nodata, /noerase, $
position=[x1,y1+(3.0-1)*b*xsize/ysize,x2,y1+(3.0+0)*b*xsize/ysize], $
xtitle=' ',   ytitle=' log(Luminosity) [1600-6000 A]', charsize=1.5, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[40,45.5],yticks=11, ystyle=1, xrange=xrange, xstyle=1,  $
xticks=nxticks, xtickv=xtickvalues, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), ytickname=['40',' ','41',' ','42',' ','43',' ','44',' ','45',' ']

; ytitle='log(Luminosity) [1600-6000 A]', 

for s=0,nSNe-1 do begin
n=s
;cgoplot, allsnflux_array[s,0,*]-SNexplosiondates[n], allsnlumlog_array[s,*], linestyle=0, thick=3,  color=altcolors[n]

cgoplot, allsnflux_array[s,0,*]-SNexplosiondates[n], allsnlumlog_array[s,*], linestyle=2, psym=-symcat(allsymnum[n]), symsize=allsizes[n], color=allcolors[n]

endfor

;;;;;;;;;;;;;;;;;;;


plot, xdata, ydata, /nodata, /noerase, $
position=[x1,y1+(2-1)*b*xsize/ysize,x2,y1+(2+0)*b*xsize/ysize], $
xtitle=' ',   ytitle='              NUV Flux Percentage         ', charsize=1.5, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[5,50.0], /ylog, yticks=5, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1)

for s=0,nSNe-1 do begin
n=s
;cgoplot, allsnflux_array[s,0,*]-SNexplosiondates[n], 100.0*allsnflux_array[s,2,*]/(allsnflux_array[s,1,*]+allsnflux_array[s,2,*]+allsnflux_array[s,3,*]), linestyle=0, thick=3,  color=altcolors[n]

cgoplot, allsnflux_array[s,0,*]-SNexplosiondates[n], 100.0*allsnflux_array[s,2,*]/(allsnflux_array[s,1,*]+allsnflux_array[s,2,*]+allsnflux_array[s,3,*]), linestyle=2, psym=-symcat(allsymnum[n]), symsize=allsizes[n], color=allcolors[n]

print, min(allsnflux_array[s,0,*]-SNexplosiondates[n])

endfor



plot, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y1+(1+0)*b*xsize/ysize], $
xtitle='Days from Explosion',   ytitle='MUV Flux Percentage ', charsize=1.5, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0.1,100.0], /ylog, yticks=3,  ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues,  ytickv=ytickvalues

for s=0,nSNe-1 do begin
n=s
;cgoplot, allsnflux_array[s,0,*]-SNexplosiondates[n], 100.0*allsnflux_array[s,1,*]/(allsnflux_array[s,1,*]+allsnflux_array[s,2,*]+allsnflux_array[s,3,*]), linestyle=0, thick=3,  color=altcolors[n]

cgoplot, allsnflux_array[s,0,*]-SNexplosiondates[n], 100.0*allsnflux_array[s,1,*]/(allsnflux_array[s,1,*]+allsnflux_array[s,2,*]+allsnflux_array[s,3,*]), linestyle=2, psym=-symcat(allsymnum[n]), symsize=allsizes[n], color=allcolors[n]
print, min(allsnflux_array[s,0,*]-SNexplosiondates[n])

endfor



;;;;;;;;;;;;;;;;


al_legend, SNlegend,   psym=allsymnum[0:nSNe-1], color=allcolors, symsize=allsizes,background_color='white', $
pos=[0.648,0.476], /norm, charsize=1.5, box=1

;

device, /close
SET_PLOT, 'X'
$open boloflux3plots_revised.eps 

print, 'start single plots'
;;;;;;;;; single plots



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  1 PLOT  ;;;;;;;;;;;;;;;;;;;;
;

; from http://www.iluvatar.org/~dwijn/idlfigures
!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
xsize = 8.8
wall = 0.03
leftmargin=0.18
bottommargin=0.14
a = xsize/8.8 - (leftmargin + wall)
b = a * 2d / (1 + sqrt(5))

ysize = (bottommargin + b + wall)*xsize
ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a

x1 = leftmargin*8.8/xsize
x2 = x1 + a*8.8/xsize
xc = x2 + wall*8.8/xsize
y1 = bottommargin*8.8/ysize
y2 = y1 + b*8.8/ysize
y3 = y2 + wall*8.8/ysize
y4 = y3 + b*8.8/ysize
yc = y4 + wall*8.8/ysize
fontsize=12

xdata=[1,2,3,4]
ydata=[1,2,3,4]

xrange=[-5,70]


yrangetop=[0.0,0.55]
yrangebottom=[0.0,0.12]

xdata=[1,2,3,4]
ydata=[2,3,4,5]

nxticks=15
nyticks=3

figurename='singlebolofluxplot.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=12, bits_per_pixel=8, /color


plot, xdata, ydata, /nodata, /noerase, $
position=[x1,y1,x2,y2], $
xtitle='Days from Explosion',   ytitle=' log(Luminosity) [1600-6000 A]', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[40,45.5],yticks=11, ystyle=1, xrange=xrange, xstyle=1,  $
xticks=nxticks, xtickv=xtickvalues, ytickv=ytickvalues, charsize=1


for s=0,nSNe-1 do begin
n=s
;cgoplot, allsnflux_array[s,0,*]-SNexplosiondates[n], allsnlumlog_array[s,*], linestyle=0, thick=3,  color=altcolors[n]

cgoplot, allsnflux_array[s,0,*]-SNexplosiondates[n], allsnlumlog_array[s,*], linestyle=2, psym=-symcat(allsymnum[n]), symsize=allsizes[n], color=allcolors[n]

endfor


al_legend, SNlist,   psym=allsymnum[0:nSNe-1], color=allcolors, symsize=allsizes, background_color='white',$
pos=[0.745,0.956], /norm, charsize=0.8, box=1


device, /close
SET_PLOT, 'X'
$open singlebolofluxplot.eps 

;;;;;;;;;;;;
;;;;;;;;

figurename='singlebolofluxnuv.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=12, bits_per_pixel=8, /color


;;;;;;;;;;;;;;;;
;for q=0,n_elements(spectralist2)-1 do begin

q=0



cgplot, xdata, ydata, /nodata, /noerase, $
position=[x1,y1,x2,y2], $
xtitle='Days from Explosion',   ytitle='NUV Flux Percentage', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[5,50.0], /ylog, yticks=5, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, ytickv=ytickvalues

print, 'first plot command for bolofluxnuv.eps'

for s=0,nSNe-1 do begin
n=s
;cgoplot, allsnflux_array[s,0,*]-SNexplosiondates[n], 100.0*allsnflux_array[s,2,*]/(allsnflux_array[s,1,*]+allsnflux_array[s,2,*]+allsnflux_array[s,3,*]), linestyle=0, thick=3,  color=altcolors[n]

cgoplot, allsnflux_array[s,0,*]-SNexplosiondates[n], 100.0*allsnflux_array[s,2,*]/(allsnflux_array[s,1,*]+allsnflux_array[s,2,*]+allsnflux_array[s,3,*]), linestyle=2, psym=-symcat(allsymnum[n]), symsize=allsizes[n], color=allcolors[n]

print, min(allsnflux_array[s,0,*]-SNexplosiondates[n])

endfor

;;;;;;;;;;;;;;;;


al_legend, SNlist,   psym=allsymnum[0:nSNe-1], color=allcolors, symsize=allsizes, background_color='white',$
pos=[0.745,0.956], /norm, charsize=0.8, box=1


device, /close
SET_PLOT, 'X'
$open singlebolofluxnuv.eps 



figurename='singlebolofluxmuv.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=12, bits_per_pixel=8, /color


;;;;;;;;;;;;;;;;
;for q=0,n_elements(spectralist2)-1 do begin

q=0


plot, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='Days from Explosion',   ytitle='MUV Flux Percentage', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0.1,100.0], /ylog, yticks=3,  ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues,  ytickv=ytickvalues, charsize=1


for s=0,nSNe-1 do begin
n=s
;cgoplot, allsnflux_array[s,0,*]-SNexplosiondates[n], 100.0*allsnflux_array[s,1,*]/(allsnflux_array[s,1,*]+allsnflux_array[s,2,*]+allsnflux_array[s,3,*]), linestyle=0, thick=3,  color=altcolors[n]

cgoplot, allsnflux_array[s,0,*]-SNexplosiondates[n], 100.0*allsnflux_array[s,1,*]/(allsnflux_array[s,1,*]+allsnflux_array[s,2,*]+allsnflux_array[s,3,*]), linestyle=2, psym=-symcat(allsymnum[n]), symsize=allsizes[n], color=allcolors[n]
print, min(allsnflux_array[s,0,*]-SNexplosiondates[n])

endfor


al_legend, SNlist,   psym=allsymnum[0:nSNe-1], color=allcolors, symsize=allsizes, background_color='white',$
pos=[0.745,0.956], /norm, charsize=0.8, box=1


device, /close
SET_PLOT, 'X'
$open singlebolofluxmuv.eps 




print, 'final stop'
stop
end
