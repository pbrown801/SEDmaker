# Imports
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import math
import scipy.interpolate
from utilities import filterlist_to_filterfiles, specarray_to_counts, isfloat, sp


def SED(spectrum, filter_list, flux_convert):
    '''
    Parameters: 
        spectrum      -- The spectrum file you wish to overlay the SED on.
        filter_list  -- The filters you wish to create the SED points for.
        flux_convert -- flux converion numbers for each filter_list in the same order.

    Use the filterlist_to_filterfiles function to get the filter files from the filters folder,
    zeropoint values, and the pivot wavelength for each filter.

    Use specarray_to_counts to get the counts for each filter. 
        - Calculate counts based on the spectrum wavelength vs flux and filter wavelength vs effectiveAreas 

    Convert the count rate to flux density using the flux conversion constant and dividing by + 10^16
    These flux densities correspond to the pivot wavelegnth of each filter.

    Plot these points on top of the spectrum. The plot is wavelength vs. Flux density.

    '''

    # Creates a list of filter files, list of zeropoint values, and list of pivot wavelengths corresponding to the filter list inputted from the parameters.
    filter_file_list, zeropointlist, pivotlist = filterlist_to_filterfiles(
        filter_list, spectrum)

    # Gets spectrum data and separates it into two lists for spectra wavelength and spectra flux.
    spectrum_data = list(map(sp, open('../spectra/'+spectrum).readlines()))
    spectrum_wave = []
    spectrum_flux = []
    for i in range(len(spectrum_data)):
        if isfloat(spectrum_data[i][0]) or isfloat(spectrum_data[i][1]):
            spectrum_wave.append(float(spectrum_data[i][0]))
            spectrum_flux.append(float(spectrum_data[i][1]))

    # Convert the lists to numpy arrays.
    spectrum_wavelength = np.array(spectrum_wave)
    spectrum_flux_dens = np.array(spectrum_flux)
    spectrum_orig = np.column_stack(
        (spectrum_wavelength[:], spectrum_flux_dens[:]))

    # Pass the spectrum information along with the filter file list to specarray_to_counts to get counts (count rates) for each filter.
    count = specarray_to_counts(spectrum_orig, filter_file_list)

    # Convert count to flux using flux conversion values
    flux = []
    for idx, i in enumerate(flux_convert):
        flux.append((count[idx]*i)/(pow(10, 16)))

    # Interpolate between the pivot wavelength and flux density
    f = scipy.interpolate.interp1d(
        pivotlist, flux, kind='linear', fill_value="extrapolate")

    # Wavelength range 1600-8000 Angstroms
    extrap = [i*10 for i in range(160, 801)]

    # Using the fill_value parameter of scipy interpolate we can extrapolate more values outside the pivot wavelength range
    flux_interp = f(extrap)
    spectrum_interp = np.column_stack((extrap, flux_interp))

    # Pass the new interpolated spectrum to specarray_to_counts
    count_interp = specarray_to_counts(spectrum_interp, filter_file_list)

    # Convert count to flux using flux conversion values
    flux_interp = []
    for idx, i in enumerate(flux_convert):
        flux_interp.append((count_interp[idx]*i)/(pow(10, 16)))

    print("FILTER LIST: ", filter_list)
    print("COUNT: ", count)
    print("INTERP COUNT: ", count_interp)
    print("FLUX: ", flux)
    print("INTERP FLUX: ", flux_interp)

    error = []
    for idx, i in enumerate(count):
        error.append((np.abs(i-count_interp[idx])/i)*100)
    print("PERCENT ERROR: ", error)

    # Ouput the pivot wavelengths values and corresponding Flux density to csv file.
    with open('../output/csv/'+spectrum[: spectrum.find(".")]+'_info.csv', mode='w') as csv_file:
        writer = csv.writer(csv_file, delimiter=",",
                            quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["Pivot Wavelength", "Flux Density"])
        for idx, i in enumerate(pivotlist):
            writer.writerow([i, flux[idx]])

    # Plotting pivot wavelength vs flux density along with the interpolated flux densities
    plt.plot(spectrum_wave, spectrum_flux)  # Spectrum data
    plt.plot(pivotlist, flux, 'r--')  # Red dashes to see the curve
    plt.plot(pivotlist, flux, 'ro')  # Red points to see the particular points
    plt.plot(pivotlist, flux_interp, 'g--')
    plt.plot(pivotlist, flux_interp, 'go')
    plt.xlim(0, pivotlist[-1]+12500)
    plt.xlabel("Wavelength (Angstroms)")
    plt.ylabel("Flux Density")
    plt.title(
        "SED for " + spectrum[: spectrum.find(".")].capitalize()+" Spectrum")
    plt.savefig('../output/plots/'+spectrum[: spectrum.find(".")]+'_SED.png')
    plt.show()


def main():
    filter_file_list = ['UVW2', 'UVM2', 'UVW1',  'U', 'B', 'V']
    flux_conversion = [6.03, 8.30, 4.02, 1.44, 1.16, 2.62]
    spectrum = 'vega.dat'
    SED(spectrum, filter_file_list, flux_conversion)


if __name__ == "__main__":
    main()
