# Imports
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import math
import scipy.interpolate
import os
from utilities import filterlist_to_filterfiles_2, specarray_to_counts, isfloat, sp


def high_res_SED(pivotlist, count, flux, filter_file_list, flux_convert):
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

    # percent error for the count rate points for the original count rates from the inputted spectrum.
    error = []
    for idx, i in enumerate(count):
        error.append((np.abs(i-count_interp[idx])/i)*100)
    print("PERCENT ERROR: ", error)

    return flux_interp

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
    # Spectrum check in spectra folder
    spectrum_dir = ''
    walk = os.walk('../spectra')
    for i in walk:
        for j in range(len(i)):
            if spectrum in i[j]:
                spectrum_dir = i[j-2]
    if spectrum_dir != '':
        spectrum_dir = (spectrum_dir.replace('../spectra', '') +
                        "\\"+spectrum).replace('\\', '/')
        # Creates a list of filter files, list of zeropoint values, and list of pivot wavelengths corresponding to the filter list inputted from the parameters.
        filter_file_list, pivotlist = filterlist_to_filterfiles_2(
            filter_list, spectrum_dir[1:])

        # Gets spectrum data and separates it into two lists for spectra wavelength and spectra flux.
        spectrum_data = list(
            map(sp, open('../spectra/'+spectrum_dir[1:]).readlines()))
        spectrum_wave = []
        spectrum_flux = []
        for i in range(len(spectrum_data)):
            if len(spectrum_data[i]) != 1:
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

        # Convert count to flux using flux conversion values from input
        flux = []
        for idx, i in enumerate(flux_convert):
            flux.append((count[idx]*i)/(pow(10, 16)))

        # Calculate avg flux +/- 50 Angstroms from pivotwavelengths of each filter
        df = pd.DataFrame({'spectrum_wave': spectrum_wave, 'spectrum_flux': spectrum_flux})
        flux_avg = []
        for p in pivotlist:
            temp=(df.loc[(df['spectrum_wave']>=p-50) & (df['spectrum_wave']<=p+50)])
            flux_avg.append(temp['spectrum_flux'].sum()/len(temp['spectrum_flux']))

        # Calculate conversion factors avg flux / count rate
        new_flux_convert = []
        for idx, f in enumerate(flux_avg):
            new_flux_convert.append((f/count[idx])*pow(10,16))
        print("ORIG CONVERSION FACTORS: ", flux_convert)
        print("NEW CONVERSION FACTORS: ",new_flux_convert)

        flux_2 = []
        for idx, i in enumerate(new_flux_convert):
            flux_2.append((count[idx]*i)/(pow(10, 16)))

        # Not changing conversion factors but using the 6 flux points to make a extrapolated spectrum and get new count rates 
        # which are then converted to flux 
        flux_interp = high_res_SED(pivotlist, count, flux, filter_file_list, flux_convert)
        # Using conversion factors created we run it through high_res SED to get new count rates. 
        flux_interp_new_convert = high_res_SED(pivotlist, count, flux_2, filter_file_list, new_flux_convert)

        # Ouput the pivot wavelengths values and corresponding Flux density to csv file.
        with open('../output/csv/'+spectrum[: spectrum.find(".")]+'_info.csv', mode='w') as csv_file:
            writer = csv.writer(csv_file, delimiter=",",
                                quotechar='"', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["Pivot Wavelength", "Flux Density"])
            for idx, i in enumerate(pivotlist):
                writer.writerow([i, flux[idx]])

        # Plotting pivot wavelength vs flux density along with the interpolated flux densities
        plt.plot(spectrum_wave, spectrum_flux)  # Spectrum data
        # Points from initial conversion of count rates to flux 
        plt.plot(pivotlist, flux, 'r--')  # Red dashes to see the curve
        # Red points to see the particular points
        plt.plot(pivotlist, flux, 'ro')
        # points from high_res_SED funcs using the spectrum extrapolated from the initial 6 points.
        plt.plot(pivotlist, flux_interp, 'g--')
        plt.plot(pivotlist, flux_interp, 'go')
        plt.plot(pivotlist, flux_interp_new_convert, 'y--')
        plt.plot(pivotlist, flux_interp_new_convert, 'yo')
        # Plot settings
        plt.xlim(0, pivotlist[-1]+12500)
        plt.xlabel("Wavelength (Angstroms)")
        plt.ylabel("Flux Density")
        plt.title(
            "SED for " + spectrum[: spectrum.find(".")].capitalize()+" Spectrum")
        plt.savefig('../output/plots/' +
                    spectrum[: spectrum.find(".")]+'_SED.png')
        plt.show()


def main():
    filter_file_list = ['UVW2', 'UVM2', 'UVW1',  'U', 'B', 'V']
    flux_conversion = [6.03, 8.30, 4.02, 1.44, 1.16, 2.62]
    spectra = ['SN2017erp_m1_UVopt.dat']
    for spectrum in spectra:
        SED(spectrum, filter_file_list, flux_conversion)


if __name__ == "__main__":
    main()
