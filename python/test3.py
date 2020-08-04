# Imports
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import math
import scipy.interpolate
import os
from utilities import filterlist_to_filterfiles_2, specarray_to_counts, isfloat, sp


def spec_dir(spectrum):
    '''
    Parameter: 
        sprectum  -- The spectrum file
    Returns the spectrum directory in the spectra folder and a boolean value to show if the file exists in the spectra folder.
    '''
    spectrum_dir = ''
    spectrum_bool = False
    walk = os.walk('../spectra')
    for i in walk:
        for j in range(len(i)):
            if spectrum in i[j]:
                spectrum_dir = i[j-2]
    if spectrum_dir != '':
        spectrum_dir = (spectrum_dir.replace('../spectra', '') +
                        "\\"+spectrum).replace('\\', '/')
        spectrum_bool = True
    return spectrum_dir, spectrum_bool


def spec_data(spectrum_dir):
    '''
    Parameter: 
        sprectum  -- The spectrum file
    Gets spectrum data and separates it into two lists: spectra wavelength and spectra flux.
    Convert the two lists into numpy arrays and make a numpy 2d array for specarray function
    Convert the two lists into dataframe 
    Returns the numpy array and dataframe
    '''
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
    spectrum_df = pd.DataFrame(
        {'spectrum_wave': spectrum_wave, 'spectrum_flux': spectrum_flux})
    return spectrum_orig, spectrum_df


def calc_convert_factors(df, pivotlist, count_orig):
    '''
    Parameters: 
        df         -- The spectrum data in a pandas dataframe
        pivotlist  -- The pivot wavelengths of the filters
        count_orig -- The original count rates from the spectrum and input conversion factors

    Calculate the average flux density +/- 50 Angstroms from each pivot wavelength.
    Divide the average flux density by original counts and multiply by + 10^16 to get calculated flux converstion factors.
    Returns the list of new flux conversion factors.
    '''
    flux_avg = []
    for p in pivotlist:
        temp = (df.loc[(df['spectrum_wave'] >= p-50)
                       & (df['spectrum_wave'] <= p+50)])
        flux_avg.append(temp['spectrum_flux'].sum()/len(temp['spectrum_flux']))

    # Calculate conversion factors avg flux / count rate
    new_flux_convert = []
    for idx, f in enumerate(flux_avg):
        new_flux_convert.append((f/count_orig[idx])*pow(10, 16))
    # print("NEW CONVERSION FACTORS: ", new_flux_convert)
    return new_flux_convert


def count_to_flux(count_orig, flux_convert):
    '''
    Parameters:
        count_orig   -- The original count rates from the spectrum and input conversion factors
        flux_convert -- The conversion factors from count rates to flux density
    Convert count to flux by multiplying the count rates by flux conversion factors.
    Returns flux
    '''
    flux = []
    for idx, i in enumerate(flux_convert):
        flux.append((count_orig[idx]*i)/(pow(10, 16)))
    return flux


def percent_error(count_orig, count_interp):
    '''
    Parameters:
        count_orig   -- The original count rates from the spectrum and input conversion factors
        count_interp -- The count rates of other conversion factors or from the high_res_sed functton 
    percent error for the count rate points for the original count rates vs interpreted counts.
    Returns percent error
    '''
    per_error = []
    for idx, i in enumerate(count_orig):
        per_error.append((np.abs(i-count_interp[idx])/i)*100)
    # print("PERCENT ERROR: ", error)
    return per_error


def high_res_SED(pivotlist, count_orig, flux, filter_file_list, flux_convert):
    '''
     Parameters: 
        pivotlist         -- The pivot wavelengths of the filters
        count_orig        -- The original count rates from the spectrum and input conversion factors
        flux              -- Flux densities converted from the count rates
        filter_file_list  -- List of filter file directories for the filters you wish to create the SED points.
        flux_convert      -- flux converion numbers for each filter_list in the same order.

    Extrapolating the flux density against the wavelength range of 1600 to 8000 Angstroms
    Using this new spectrum to generate new count rates using specarray_to_counts function
    Convert the count rates to flux using the inputed flux_convert values.
    Calculate percent error using percent_error function.
    Returns the new extrapolated flux and error found.
    '''
    # Interpolate between the pivot wavelength and flux density
    f = scipy.interpolate.interp1d(
        pivotlist, flux, kind='linear', fill_value="extrapolate")

    # Wavelength range 1600-8000 Angstroms
    extrap = [i*10 for i in range(160, 801)]

    # Using the fill_value parameter of scipy interpolate we can extrapolate more values outside the pivot wavelength range
    flux_extrap = f(extrap)
    spectrum_interp = np.column_stack((extrap, flux_extrap))

    # Pass the new interpolated spectrum to specarray_to_counts
    count_extrap = specarray_to_counts(spectrum_interp, filter_file_list)

    # Convert count to flux using flux conversion values
    flux_extrap = count_to_flux(count_extrap, flux_convert)

    # percent error for the count rate points for the original count rates from the inputted spectrum.
    error = percent_error(count_orig, count_extrap)
    return flux_extrap, error


def SED_plot(spectrum, spectrum_df, pivotlist, flux_convert_name):
    '''
    Parameters: 
        spectrum     -- The spectrum file you wish to overlay the SED on.
        spectrum_df  -- The spectrum data (wavelength and flux densities)
        pivotlist    -- The pivot wavelengths of the filters

    Set up the figure and pyplot plotting mechanics.
    Set up the plot settings.
    Plot the spectrum itself.
    Returns the ax plot object
    '''
    # plt.style.use('classic')

    fig, ax = plt.subplots()
    # Plot spectrum
    ax.plot(spectrum_df['spectrum_wave'], spectrum_df['spectrum_flux'], label='Spectrum')
    # Plot Settings
    ax.set_xlabel("Wavelength (Angstroms)")
    ax.set_ylabel("Flux Density")
    ax.set_title(
        "SED for " + spectrum[: spectrum.find(".")].capitalize()+"Spectrum w/ " +flux_convert_name +" conversion input")
    ax.axis = ('equal')
    return ax


def SED(spectrum, filter_list, flux_convert, flux_conversion_static):
    '''
    Parameters: 
        spectrum     -- The spectrum file you wish to overlay the SED on.
        filter_list  -- The filters you wish to create the SED points for.
        flux_convert -- flux converion numbers for each filter_list in the same order.
        flux_conversion_static -- flux conversion factors used with every spectrum that is inputted. 
                                  If vega or Pickles is the main flux_convert then it will not be created

    Use the spec_dir function to check if spectrum file exists.

    Use the filterlist_to_filterfiles function to get the filter files from the filters folder,
    zeropoint values, and the pivot wavelength for each filter.

    Use specarray_to_counts to get the counts for each filter. 
        - Calculate counts based on the spectrum wavelength vs flux and filter wavelength vs effectiveAreas 

    Convert the count rate to flux density using the flux conversion constant and dividing by + 10^16
    These flux densities correspond to the pivot wavelegnth of each filter.

    Plot these points on top of the spectrum. The plot is wavelength vs. Flux density.
    '''
    # Need to extract the conversion input name and values
    for key in flux_convert:
        flux_convert_name = key
        flux_convert=flux_convert[flux_convert_name]

    # check if spectrum file is in spectra directory
    spectrum_dir, spectrum_bool = spec_dir(spectrum)
    
    if spectrum_bool:

        # Creates a list of filter files, list of zeropoint values, and list of pivot wavelengths corresponding to the filter list inputted from the parameters.
        filter_file_list, pivotlist = filterlist_to_filterfiles_2(
            filter_list, spectrum_dir[1:])
        # Gets the spectrum data from the spectrum file
        spectrum_orig, spectrum_df = spec_data(spectrum_dir)

        # Counts rate from the spectrum data
        count_orig = specarray_to_counts(spectrum_orig, filter_file_list)

        # Calculating new conversion factors based on avg flux (+-50 Angstroms) / count_orig
        new_flux_convert = calc_convert_factors(
            spectrum_df, pivotlist, count_orig)

        # Convert the count rates to flux densities using various conversions
        flux_comparison = count_to_flux(count_orig, flux_convert)

        if flux_convert != flux_conversion_static[0]:
            flux_vega = count_to_flux(count_orig, flux_conversion_static[0])
            flux_extrap_vega, error_extrap_vega = high_res_SED(
                pivotlist, count_orig, flux_vega, filter_file_list, flux_conversion_static[0])
        else:
            flux_extrap_vega, error_extrap_vega = high_res_SED(
                pivotlist, count_orig, flux_comparison, filter_file_list, flux_conversion_static[0])
        
        if flux_convert != flux_conversion_static[1]:
            flux_Pickles = count_to_flux(count_orig, flux_conversion_static[1])
        
        flux_calc_conversion = count_to_flux(count_orig, new_flux_convert)

        # Plotting - Uses the ax plot object to plot the flux densities found from using different methods
        ax = SED_plot(spectrum, spectrum_df, pivotlist, flux_convert_name)
        ax.plot(pivotlist, flux_comparison, 'ro--', label='comparison')
        if flux_convert != flux_conversion_static[0]:
            ax.plot(pivotlist, flux_vega, 'go--', label='vega')
        ax.plot(pivotlist, flux_extrap_vega, 'yo--', label='vega+interation')
        if flux_convert != flux_conversion_static[1]:
            ax.plot(pivotlist, flux_Pickles, 'mo--', label='pickles')
        ax.plot(pivotlist, flux_calc_conversion, 'ko--', label='calc convert')
        ax.legend()
        plt.savefig('../output/plots/' +
                    spectrum[: spectrum.find(".")]+'_SED_with'+flux_convert_name+'_conversion.png')
        plt.show()

    else:
        print("SPECTRUM FILE DOES NOT EXIST.")


def main():
    flux_conversion_dict = {
        'vega':  [6.03, 8.30, 4.02, 1.44, 1.16, 2.62], 'GRBs': [5.98, 8.45, 4.21, 1.63, 1.47, 2.61], 'Pickles': [5.77, 7.47, 4.06, 1.53, 1.31, 2.61],
        'AB': [6.23, 8.49, 4.63, 1.66, 1.48, 2.61], 'ST': [6.03, 8.30, 4.02, 1.44, 1.16, 2.62], '3000 K': [0.03, 2.53, 0.57, 1.49, 1.31, 2.53],
        '10000 K': [5.76, 8.28, 4.55, 1.64, 1.49, 2.62], '30000 K': [6.06, 8.42, 4.02, 1.57, 1.48, 2.62], 'a0i': [6.24, 7.53, 4.01, 1.37, 1.34, 2.63],
        'a0iii': [5.86, 7.56, 4.42, 1.45, 1.14, 2.62], 'a0v': [6.05, 7.94, 4.15, 1.44, 1.19, 2.62], 'g0v': [0.09, 6.72, 2.89, 1.61, 1.33, 2.60],
        'o9v': [6.08, 7.96, 4.29, 1.58, 1.40, 2.58], 'g1050 04': [6.97, 8.87, 3.55, 1.37, 1.67, 2.51], 'ic3639': [6.09, 6.97, 3.82, 1.47, 1.52, 2.53],
        'mrk477': [6.83, 8.98, 4.32, 1.32, 1.22, 1.46], 'ngc6221': [1.61, 5.23, 3.83, 1.59, 1.40, 2.69], 'ngc7496': [7.42, 7.90, 4.17, 1.49, 1.40, 2.66],
        'IC 4051': [2.50, 8.93, 2.94, 1.65, 1.30, 2.57], 'IC 5298': [5.75, 8.07, 3.99, 1.53, 1.36, 2.57], 'II Zw 096': [6.27, 8.54, 4.24, 1.45, 1.58, 2.64],
        'NGC 0520': [5.52, 8.44, 3.84, 1.61, 1.31, 2.62], 'NGC 0584': [2.30, 8.92, 2.96, 1.56, 1.37, 2.63], 'Hsiao 0': [3.67, 6.81, 1.91, 2.15, 1.14, 2.28],
        'Hsiao 15': [2.52, 6.11, 1.24, 1.65, 1.28, 3.13], 'SN05cs+3': [2.81, 9.11, 4.61, 1.58, 1.55, 2.57], 'SN05cs+17': [1.11, 1.40, 2.70, 3.13, 1.26, 2.78],
        'Ia SN2011fe': [3.20, 2.95, 0.93, 2.01, 1.14, 2.42], 'Ia SN1992A': [3.82, 3.32, 0.85, 1.75, 1.08, 2.48], 'Ic SN1994I': [3.40, 7.53, 1.50, 1.91, 1.16, 2.96],
        'IIP SN1999em': [5.59, 6.45, 3.75, 1.69, 1.68, 1.33]}
    filter_file_list = ['UVW2', 'UVM2', 'UVW1',  'U', 'B', 'V']
    flux_convert = {'vega': flux_conversion_dict['vega']}
    flux_conversion_static = [
        flux_conversion_dict['vega'], flux_conversion_dict['Pickles']]
    spectra = ['SN2017erp_m1_UVopt.dat']
    for spectrum in spectra:
        SED(spectrum, filter_file_list, flux_convert, flux_conversion_static)


if __name__ == "__main__":
    main()
