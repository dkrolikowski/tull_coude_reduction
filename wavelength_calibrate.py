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
import os
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

def find_arc_lamp_line_pixel_centers( flux, config ):
    """ Function to find lines in an input arc lamp spectrum using scipy.signal's find_peaks algorithm, with options set in the config file.

    Parameters
    ----------
    flux : array
        The flux spectrum array.
    config : dict
        The overall config file defined using YAML with all parameters for running the reduction and analysis pipeline.

    Returns
    -------
    peak_pixels_initial_fit : array
        The array of pixel centroids found as lines in the input flux spectrum.
    """
    
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
        try:
            line_fit, _ = optimize.curve_fit( tull_coude_utils.gaussian_1d, fit_range, flux[fit_range], p0 = p_guess )
        except:
            continue
        
        # Output the new fit centroid
        peak_pixels_initial_fit[i_peak] = line_fit[1]
        
    # Get rid of any nans!
    peak_pixels_initial_fit = peak_pixels_initial_fit[np.isfinite( peak_pixels_initial_fit )]
        
    return peak_pixels_initial_fit

def fit_wavelength_solution( pixel_centroids, prelim_wavelengths, line_list_wavelengths, config ):
    
    wavsol = np.zeros( len(spec) ) # Initialize the wavelength solution for this order

    # Find peaks in the arc spectrum!
    pixcent, wavcent = Find_Peaks( wav, spec, specsig, peaksnr = snr, minsep = minsep )

    # Initialize kept/rejected peaks, and match to nearest in the THAR catalog
    keeps = { 'pix': np.array([]), 'wav': np.array([]), 'line': np.array([]) }
    rejs  = { 'pix': np.array([]), 'wav': np.array([]), 'line': np.array([]) }

    for i in range( len( wavcent ) ): # Find matches to peaks in the THAR catalogue
        dists    = np.absolute( THARcat - wavcent[i] )
        mindist  = np.argmin( dists )

        if dists[mindist] <= 1.0: # If the wavelengths of the lines are close keep them!
            keeps['pix']  = np.append( keeps['pix'], pixcent[i] )
            keeps['wav']  = np.append( keeps['wav'], wavcent[i] )
            keeps['line'] = np.append( keeps['line'], THARcat[mindist] )

    # Now actually do the fit!
    dofit  = True
    ploti  = 1
    cutoff = 4.0 / 0.67449 # Corrects MAD to become sigma (x4)

    while dofit:

        # Polynomial fit to wavelength solution
        res = np.polyfit( keeps['pix'], keeps['line'], Conf.WavPolyOrd, full = True )

        wavparams  = res[0]
        fitresult  = res[1:]
        ptsfromfit = np.polyval( wavparams, keeps['pix'] )

        # Calculate residuals in wavelength and velocity
        wavresids  = ptsfromfit - keeps['line']
        velresids  = wavresids / keeps['line'] * 3e5
        resids     = { 'wav': wavresids, 'vel': velresids }
        medabsdev  = np.median( np.abs( np.abs( resids['wav'] ) - np.median( np.abs( resids['wav'] ) ) ) )

        # Determine which lines are outliers in wavelength and velocity residuals
        velcut = np.sum( np.abs(resids['vel']) >= 5.0 )

        # Reject points that are further away than the sigma cutoff
        torej  = ( np.abs( resids['wav'] ) >= cutoff * medabsdev ) # | ( keeps['pix'] < 512 ) | ( keeps['pix'] > 1536 )

        tokeep = np.logical_not( torej )
        numrej = np.sum( torej )

        if velcut > 0 and numrej != len(torej): # If there are points that are outliers!
            if numrej > 0: # If there are points to reject based on wavelength residual

                plotname = path + '/resids_round_' + str(ploti) + '.pdf'
                Plot_WavSol_Resids( resids, keeps['line'], cutoff, plotname, tokeep = tokeep, toreject = torej )

                pltwav = np.polyval( wavparams, np.arange( len( spec ) ) )
                plt.clf()
                plt.plot( pltwav, np.log10(spec), 'k-', lw = 1 )
                plt.plot( THAR['wav'], THAR['logspec'], 'r-', lw = 1 )
                plt.xlim( pltwav[0], pltwav[-1] )
                for peak in keeps['line']:
                    plt.axvline( x = peak, color = 'b', ls = ':', lw = 1 )
                plt.savefig(path + '/rejplots/fullspec_' + str(ploti) + '.pdf'); plt.clf()

                rejs['pix']  = keeps['pix'][torej]
                rejs['wav']  = keeps['wav'][torej]
                rejs['line'] = keeps['line'][torej]

                keeps['pix']  = keeps['pix'][tokeep]
                keeps['wav']  = keeps['wav'][tokeep]
                keeps['line'] = keeps['line'][tokeep]

                ploti += 1

            # If there aren't points to reject via wavelength residual, but still velocity outliers, reduce the cutoff limit by 1 sigma
            elif numrej == 0 and velcut > 0:
                cutoff = cutoff - 1.0 / 0.67449

            else: # Honestly a little unsure what this is doing....
                if Conf.verbose: print( 'There is something seriously wrong.\n' )
                if Conf.verbose: print( 'There are points > 0.2 km/s, but none are found to be rejected. FIX' )
                flag = True
                if plots:
                    plotname = path + '/resids_round_' + str(ploti) + '_flag.pdf'
                    Plot_WavSol_Resids( resids, keeps['line'], cutoff, plotname )
                break

        # If it wants to reject all the points! That's bad!
        elif numrej == len(torej):
            if plots:
                plotname = path + '/resids_round_' + str(ploti) + '.pdf'
                Plot_WavSol_Resids( resids, keeps['line'], cutoff, plotname, toreject = torej )
            flag = True
            dofit = False

        # If it all works out fine!
        else:
            if plots:
                plotname = path + '/resids_round_' + str(ploti) + '.pdf'
                Plot_WavSol_Resids( resids, keeps['line'], cutoff, plotname, tokeep = tokeep )
            flag = False
            dofit = False

    # Basically if there aren't enough peaks and the fit isn't well constrained
    if fitresult[0].size == 0:
        flag = True

    # Full wavelength solution for the order
    wavsol = np.polyval( wavparams, np.arange( len(spec) ) )

    
    
    return wavelength_solution_coeffs

##### Main wrapper script for wavelength calibration

def wavelength_solution_and_calibrate( arc_file_indices, header_df, config ):
    
    ##### First make the wavelength solution!
    
    ### Go through each of the input file indices for the arc lamps to use
    for i_file in arc_file_indices:
        
        file_in = fits.open( os.path.join( config['paths']['reduction_dir'], 'spectrum_files', 'tullcoude_{}_spectrum.fits'.format( header_df['file_token'].values[i_file] ) ) )
        
        # Go order by order and get wavelength solution
        for order in range( config['trace']['number_of_orders'] ):
            
            # Get the peak centroids
            lamp_line_pixel_centroids = find_arc_lamp_line_pixel_centers( file_in[1].data[order], config )
            
            # Fit the wavelength solution
            order_wavelength_solution = fit_wavelength_solution( lamp_line_pixel_centroids, config )
            
            # Make any plots!

    ### Then apply the wavelength solution to all of the extracted spectra
    
    return None