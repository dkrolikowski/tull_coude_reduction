""" Functions to wavelength calibrate the spectra. Finds lines in arc lamp spectra (ThAr), matches them to a line list based on a preliminary input wavelength solution, and then iteratively fits a polynomial solution.

Created by DMK on 11/11/2023
Last updated by DMK on 11/11/2023
"""

##### Imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from astropy.io import fits
from scipy import optimize, signal, stats

import glob
import tqdm

import tull_coude_utils

##### Functions

def get_flux_mad_from_spectral_chunks( flux, chunk_size = 50 ):
        
    # The number of chunks from the input chunk size in pixels
    number_of_chunks = flux.size // chunk_size
    
    # The array to hold the MAD for each of the spectral chunks
    mad_arr = np.full( number_of_chunks, np.nan )
    
    # Go through each of the chunks
    for i in range( number_of_chunks ):
    
        # The indices of the spectrum within this chunk
        ind_low, ind_high = chunk_size * i, chunk_size * ( i + 1 )
    
        # Calculate the scaled MAD
        mad_arr[i] = stats.median_abs_deviation( flux[ind_low:ind_high], nan_policy = 'omit', scale = 'normal' )
    
    # Adopt the median of the chunk MAD array
    median_chunk_mad = np.nanmedian( mad_arr )

    return median_chunk_mad

def find_arc_lamp_lines( flux, config ):
    
    ### First prepare from input flux spectrum: normalize it and get noise estimate
    
    # Normalize the flux by its median
    flux = flux / np.nanmedian( flux )
    
    # Now get an estimate of the noise for peak finding threshold. Two current options: one is just from the MAD of the flux spectrum, the other divides the spectrum into chunks and takes the median of each chunk's MAD.
    if config['wavecal']['peak_threshold_mad_method'] == 'full_spectrum':
        
        mad_peak_threshold = stats.median_abs_deviation( flux, nan_policy = 'omit', scale = 'normal' )

    elif config['wavecal']['peak_threshold_mad_method'] == 'chunk_spectrum':
        
        mad_peak_threshold = get_flux_mad_from_spectral_chunks( flux )
        
    ### Now find the peaks using scipy signal's find peak algorithm, with the prominence threshold based on the "noise" estimate of the spectrum
    
    peak_pixels_initial, peak_properties = signal.find_peaks( flux, distance = config['wavecal']['lamp_line_min_separation_pix'], width = config['wavecal']['lamp_line_pix_width_limits'], prominence = config['wavecal']['lamp_line_peak_threshold_sigma'] * mad_peak_threshold )
    
    # Now go through each of the lines found and fit with a Gaussian to get the decimal pixel centroid
    
    peak_pixels_initial_fit = np.full( peak_pixels_initial.size, np.nan )
    
    for i_peak, peak in enumerate( peak_pixels_initial ):
        
        # Get the fit location: within +/- the config 'lamp_line_min_separation_pix' of the peak centroid
        fit_range = np.arange( peak - config['wavecal']['lamp_line_min_separation_pix'], peak + config['wavecal']['lamp_line_min_separation_pix'] + 1 )
        # Get rid of any nans!
        fit_range = np.intersect1d( fit_range, np.where( np.isfinite( flux ) )[0] )

        # The initial guess for 1D gaussian line parameters: the flux at the found centroid, the found int centroid, the minimum of the config 'lamp_line_pix_width_limits', and 0
        p_guess = [ flux[peak], peak, min( config['wavecal']['lamp_line_pix_width_limits'] ), 0 ]
        
        # Fit!
        line_fit, _ = optimize.curve_fit( tull_coude_utils.gaussian_1d, fit_range, flux[fit_range], p0 = p_guess )
        
        # Output the new fit centroid
        peak_pixels_initial_fit[i_peak] = line_fit[1]
        
    return None

##### Main wrapper script for wavelength calibration