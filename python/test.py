# IMPORTS
import matplotlib.pyplot as plt
import seaborn as sns
import time
import pandas as pd
import numpy as np
import math
from utilities import filterlist_to_filterfiles, specarray_to_counts, isfloat, sp

# filter_file_list = ['UVW2', 'UVM2', 'UVW1',  'U', 'B', 'V', 'R', 'I']
filter_file_list = ['UVW2', 'UVM2', 'UVW1',  'U', 'B', 'V']
flux_conversion = [6.03, 8.30, 4.02, 1.44, 1.16, 2.62]
filter_file_list,zeropointlist,pivotlist = filterlist_to_filterfiles(filter_file_list, 'vega.dat')
spectra_data = list(map(sp, open('../spectra/vega.dat').readlines()))

spectra_wave = []
spectra_flux = []
for i in range(len(spectra_data)):
    if isfloat(spectra_data[i][0]) or isfloat(spectra_data[i][1]):
        spectra_wave.append(float(spectra_data[i][0]))
        spectra_flux.append(float(spectra_data[i][1]))

spectra_wavelength = np.array(spectra_wave)
spectra_flux_dens = np.array(spectra_flux)
spectrum = np.column_stack((spectra_wavelength[:],spectra_flux_dens[:]))
count = specarray_to_counts(spectrum, filter_file_list)
flux = []
for idx,i in enumerate(flux_conversion):
    flux.append((count[idx]*i)/(pow(10, 16)))
print("COUNT: ",count)
print("Flux: ",flux)
print("Pivot: ",pivotlist)

fig=plt.figure()
ax=plt.axes()
ax.plot(spectra_wave, spectra_flux)
ax.plot(pivotlist, flux, marker='.')
plt.xlim(1000, 20000)
plt.xlabel("Wavelength")
plt.ylabel("Flux Density")
plt.title("SED for Vega spectrum")
# plt.show()

# trying to find bins to automate the xlim of plot 
div = [i*100 for i in range(1, 151)]
for d in div:
    # print(d)
    bins = math.ceil((3000000-900.5)/d)-1
    if bins < 300:
        print(d)
        print(bins)
        break
wavelength_bins = np.empty(bins)
for i in range(bins):
    wavelength_bins[i]=900.5+d*i
print(wavelength_bins)