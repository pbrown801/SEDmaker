# Imports
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
from utilities import filterlist_to_filterfiles, specarray_to_counts, isfloat, sp


def SED(spectra, filter_list, flux_convert):
    '''
    Parameters: 
        spectra      -- The spectrum file you wish to overlay the SED on.
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

    #Creates a list of filter files, list of zeropoint values, and list of pivot wavelengths corresponding to the filter list inputted from the parameters. 
    filter_file_list, zeropointlist, pivotlist = filterlist_to_filterfiles(filter_list, spectra)

    #Gets spectrum data and separates it into two lists for spectra wavelength and spectra flux.
    spectra_data = list(map(sp, open('../spectra/'+spectra).readlines()))
    spectra_wave = []
    spectra_flux = []
    for i in range(len(spectra_data)):
        if isfloat(spectra_data[i][0]) or isfloat(spectra_data[i][1]):
            spectra_wave.append(float(spectra_data[i][0]))
            spectra_flux.append(float(spectra_data[i][1]))

    # Convert the lists to numpy arrays.
    spectra_wavelength = np.array(spectra_wave)
    spectra_flux_dens = np.array(spectra_flux)
    spectrum = np.column_stack((spectra_wavelength[:], spectra_flux_dens[:]))

    # Pass the spectrum information along with the filter file list to specarray_to_counts to get counts (count rates) for each filter.
    count = specarray_to_counts(spectrum, filter_file_list)

    # Convert count to flux using flux conversion values
    flux = []
    for idx, i in enumerate(flux_convert):
        flux.append((count[idx]*i)/(pow(10, 16)))

    # Ouput the pivot wavelengths values and corresponding Flux density to csv file.
    with open('../output/csv/'+spectra[: spectra.find(".")]+'_info.csv', mode='w') as csv_file:
        writer = csv.writer(csv_file, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["Pivot Wavelength","Flux Density"])
        for idx,i in enumerate(pivotlist):
            writer.writerow([i, flux[idx]])

    # Plotting pivot wavelength vs flux density
    plt.plot(spectra_wave, spectra_flux) #Spectrum data
    plt.plot(pivotlist, flux, 'r--') # Red dashes to see the curve
    plt.plot(pivotlist, flux, 'ro') # Red points to see the particular points
    plt.xlim(1000, 20000)
    plt.xlabel("Wavelength")
    plt.ylabel("Flux Density")
    plt.title("SED for Vega spectrum")
    plt.savefig('../output/plots/'+spectra+'_SED.png')
    # plt.show()


def main():
    filter_file_list = ['UVW2', 'UVM2', 'UVW1',  'U', 'B', 'V']
    flux_conversion = [6.03, 8.30, 4.02, 1.44, 1.16, 2.62]
    spectra = 'vega.dat'
    SED(spectra, filter_file_list, flux_conversion)


if __name__ == "__main__":
    main()
