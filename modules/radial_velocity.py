""" Functions to measure radial velocities from extracted echelle spectra.

Created by DMK on 11/20/2023

Last updated by DMK on 11/22/2023
"""

##### Imports
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from astropy import constants
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from datetime import datetime
from scipy.optimize import curve_fit
from scipy.stats import median_abs_deviation

import barycorrpy
import os
import tqdm
import saphires

##### Functions

def make_saphires_ls_file( orders_to_use, data_wavelength, ls_file_name ):
    """ Function to generate the input ls wavelength file for sahpires broadening function computation.

    Parameters
    ----------
    orders_to_use : pandas DataFrame
        The data from the file pre-defining which orders and ranges to use for the BF computation. Has columns 'order', 'wave_start', 'wave_end', and 'flag'.
    data_wavelength : array
        The data wavelength array.
    ls_file_name : str
        Name of the ls file.

    Returns
    -------
    None.
    """
    
    # Open the temporary wavelength ls file for the BF run
    output_ls_file = open( ls_file_name, 'w' )
    
    # Go through each of the orders listed in the input file
    for i_order, order in enumerate( orders_to_use['order'].values ):
        
        # Set the starting and ending wavelength, with padding on either side
        wave_start, wave_end = data_wavelength[order].min() - 10, data_wavelength[order].max() + 10
        
        # Check for flags -- set in the input file to avoid tellurics or strong stellar features (e.g. H-alpha)
        if orders_to_use['flag'].values[i_order] == 'start':
            wave_start = orders_to_use['wave_start'].values[i_order]
        elif orders_to_use['flag'].values[i_order] == 'end':
            wave_end = orders_to_use['wave_end'].values[i_order]
        
        # Write the order's info line to the output temporary ls file
        output_ls_file.write( '{} {:.2f}-{:.2f}\n'.format( order, wave_start, wave_end ) )
        
    # Close the ls file
    output_ls_file.close()
        
    return None

def bootstrap_sample_bf_rvs( bf_tar, bf_tar_spec, num_samples ):
    """ Perform bootstrap sampling of combining order by order broadening functions to get a set of RV samples for measuring the star's RV value and error.
    This is done by combining different sets of orders with replacement.

    Parameters
    ----------
    bf_tar : array
        The array of keys for the saphires tar_spec dictionary. Each key is an order.
    bf_tar_spec : dict
        The saphires-defined dictionary for each order with broadening function/spectrum information. The result of various saphires functions.
    num_samples : int
        The number of bootstrap samples to generate.

    Returns
    -------
    rv_samples : array
        The velocity centroid of the Gaussian fit to each sampled combined BF.
    """

    # Set up empty array for holding the bootstrap RV samples
    rv_samples = np.full( num_samples, np.nan )

    # Loop through number of samples
    for i_sample in range( num_samples ):

        # Get a set of random indices between 0 and number of orders with BFs, defines the sampled-with-replacement set of orders for combining
        orders_to_sample = np.random.randint( 0, len( bf_tar ), len( bf_tar ) )

        # Weight combine the broadening functions for the sampled orders
        v_comb, bf_comb, _, _ = saphires.bf.weight_combine( bf_tar[orders_to_sample], bf_tar_spec, std_perc = 0.1 )
        
        # Fit with a Gaussian! Set initial guesses
        p0 = [ np.max( bf_comb ), v_comb[np.argmax( bf_comb )], 1, np.median( bf_comb ) ]

        try: # Uppoed the maxfev because it was failing out
            fit, errs = curve_fit( saphires.utils.gaussian_off, v_comb, bf_comb, p0 = p0, maxfev = 50000 )
        except RuntimeError:
            continue

        # Add in an rv shift if has been applied to the BF computation (here it should be zero but is included for future-proofing)
        rv_samples[i_sample] = fit[1] + bf_tar_spec[bf_tar[0]]['rv_shift']

    return rv_samples

def make_rv_compiled_excel( file_indices, output_file_name, header_df, config ):
    """ Function to read in science frames with RV measurements and output them as a CSV of the night for easy access.

    Parameters
    ----------
    file_indices : array
        The file indices (in the header information file) of science observations with extracted spectra to measure RVs for.
    output_file_name : str
        File name for the output RV compilation CSV.
    header_df : pandas DataFrame
        The compiled information from the file headers.
    config : dict
        The overall config file defined using YAML with all parameters for running the reduction and analysis pipeline.

    Returns
    -------
    None.
    """
    
    # Make output dictionary with RV and RV error arrays
    out_dict = { 'file_token': header_df['file_token'].values[file_indices], 'object': header_df['object'].values[file_indices], 'rv': np.full( file_indices.size, np.nan ), 'rv_error': np.full( file_indices.size, np.nan ) }
    
    # Go through each of the science frames
    for i_file, file_index in enumerate( file_indices ):
        
        # Read in the spectra file
        file_name = os.path.join( config['paths']['reduction_dir'], 'spectrum_files', 'tullcoude_{}_spectrum.fits'.format( header_df['file_token'].values[file_index] ) )
        file_in   = fits.open( file_name )
        
        # Make sure that the radial velocity extension exists!
        if 'radial velocity' in file_in:
            
            out_dict['rv'][i_file]       = file_in[0].header['RVBF']
            out_dict['rv_error'][i_file] = file_in[0].header['ERVBF']
            
    # Output as an excel file
    pd.DataFrame( out_dict ).to_csv( output_file_name, index = False )
    
    return None

### Plotting functions

def plot_bootstrap_rv_result( bf_tar, bf_tar_spec, bc_vel, rv_samples, file_name ):
    """ Function to make two-panel plot for the RV bootstrap result. One panel shows the order by order and combined BFs, and one panel shows the histogram of the bootstrap RV samples.

    Parameters
    ----------
    bf_tar : array
        The array of keys for the saphires tar_spec dictionary. Each key is an order.
    bf_tar_spec : dict
        The saphires-defined dictionary for each order with broadening function/spectrum information. The result of various saphires functions.
    bc_vel : float
        The barycentric velocity correction (in km/s).
    rv_samples : array
        The bootstrap RV samples.
    file_name : str
        File name to save the plot to.

    Returns
    -------
    None.
    """
    
    # First calculate the RV value and error
    rv_value = np.nanmedian( rv_samples )
    rv_error = median_abs_deviation( rv_samples, scale = 'normal', nan_policy = 'omit' )
    
    ##### Make the plot
        
    # Make the figure -- two panels
    fig, axs = plt.subplots( 1, 2, num = 1, clear = True )
    fig.set_size_inches( 15, 7 )
    
    ### Left panel: BF and Gaussian fit
    
    # Combine the BFs
    velocity_arr, combined_bf, _, _ = saphires.bf.weight_combine( bf_tar, bf_tar_spec, std_perc = 0.1 )

    # Plot the combined BF, with the barycentric correction added in
    axs[0].plot( velocity_arr + bc_vel, combined_bf, c = '#323232', lw = 1.5, label = None, zorder = 3 )
    
    # Plot the Gaussian fit!
    fit, errs = curve_fit( saphires.utils.gaussian_off, velocity_arr, combined_bf, p0 = [ np.max( combined_bf ), velocity_arr[np.argmax( combined_bf )], 1, np.median( combined_bf ) ], maxfev = 50000 )
    fit_rv    = fit[1] + bc_vel + fit[1] * bc_vel / constants.c.to('km/s').value
    
    axs[0].plot( velocity_arr + bc_vel, saphires.utils.gaussian_off( velocity_arr, *fit ), '--', c = '#bf3465', lw = 1, label = 'Fit RV: {:.3f}'.format( fit_rv ), zorder = 3 )
    
    # Plot vertical lines for the adopted RV value and 3 sigma uncertainty range
    axs[0].axvline( x = rv_value, c = '#dfa5e5', label = 'Bootstrap RV: {:.3f} +/- {:.3f}'.format( rv_value, rv_error ), lw = 1.5, zorder = 2 )
    axs[0].axvspan( rv_value - 3 * rv_error, rv_value + 3 * rv_error, color = '#dfa5e5', alpha = 0.1, zorder = 1 )
    
    # Also plot a vertical line at the barycentric correction velocity
    axs[0].axvline( x = bc_vel, c = '#874310', lw = 0.75, zorder = 2, label = 'Barycentric Velocity Correction: {:.2f} km/s'.format( bc_vel ) )

    axs[0].set_xlabel( 'Velocity (km/s)' )
    axs[0].set_ylabel( 'Brodening Function' )
    axs[0].legend( fontsize = 'x-small' )
    
    ### Right panel: RV bootstrap samples
    
    # Histogram of the bootstrap RV samples
    axs[1].hist( rv_samples, bins = 50, histtype = 'step', color = '#323232', lw = 1.5, range = ( rv_value - 10 * rv_error, rv_value + 10 * rv_error ), zorder = 3 )
    
    # Plot vertical lines for the adopted RV value and 1 sigma uncertainty range
    axs[1].axvline( x = rv_value, c = '#dfa5e5', label = 'Bootstrap RV: {:.3f} +/- {:.3f}'.format( rv_value, rv_error ), lw = 1.5, zorder = 2 )
    axs[1].axvspan( rv_value - rv_error, rv_value + rv_error, color = '#dfa5e5', alpha = 0.1, zorder = 1 )

    axs[1].set_xlabel( 'Bootstrap Sample RV (km/s)' )
    axs[1].set_ylabel( 'Histogram Counts' )
    axs[1].legend( fontsize = 'small' )

    plt.savefig( file_name, bbox_inches = 'tight', pad_inches = 0.05 )
    plt.close()
    
    return None

##### Main script

def measure_radial_velocity( file_indices, header_df, config ):
    """ Main script to run for measuring radial velocities using broadening functions.

    Parameters
    ----------
    file_indices : array
        The file indices (in the header information file) of science observations with extracted spectra to measure RVs for.
    header_df : pandas DataFrame
        The compiled information from the file headers.
    config : dict
        The overall config file defined using YAML with all parameters for running the reduction and analysis pipeline.

    Returns
    -------
    None.
    """
    
    ### Loop through each of the files to measure RV for
    for i_file in tqdm.tqdm( file_indices, desc = 'Measuring radial velocities', colour = 'green' ):
        
        ### Data preparation -- read in, normalize, get barycentric velocity correction
        
        # Read in the spectra file
        file_name = os.path.join( config['paths']['reduction_dir'], 'spectrum_files', 'tullcoude_{}_spectrum.fits'.format( header_df['file_token'].values[i_file] ) )
        file_in   = fits.open( file_name )
                
        # Continuum normalize the spectrum
        flux_cont_norm = file_in['extracted flux'].data / file_in['continuum'].data

        ## BC velocity

        # Get the data midtime. If flux weighted midtime is in the header use that, otherwise the header OBSJD + half exposure time
        if 'emfwmtim' in file_in[0].header:
            obs_mid_jd = Time( file_in[0].header['DATE-OBS'] + 'T' + file_in[0].header['emfwmtim'], format = 'isot' ).jd
        else:
            obs_mid_jd = header_df['obs_jd'].values[i_file] + 0.5 * header_df['exp_time'].values[i_file] / 60. / 60. / 24.
                        
        # Get RA/Dec from header into degrees
        sky_coord = SkyCoord( file_in[0].header['ra'], file_in[0].header['dec'], unit = ( u.hourangle, u.deg ) )
                
        # Get barycentric velocity correction from barycorrpy
        star_info    = { 'ra': sky_coord.ra.value, 'dec': sky_coord.dec.value, 'epoch': 2451545.0 }
        barycorr_vel = barycorrpy.get_BC_vel( obs_mid_jd, obsname = 'McDonald', leap_update = False, **star_info )[0][0] / 1000.0
        
        ### Saphires broadening function prep -- create ls file, make data and template structures
        
        # Read in the pre-defined orders to use file
        orders_to_use = pd.read_csv( os.path.join( config['paths']['code_dir'], 'data', config['radial_velocity']['orders_to_use_file_name'] ) )
        
        # Make the ls file
        temp_ls_file_name = 'temp_{}.ls'.format( header_df['file_token'].values[i_file] )
        make_saphires_ls_file( orders_to_use, file_in['wavelength'].data, temp_ls_file_name )
        
        # Make the saphires observed data structures. Don't combine all orders (unnecessary and makes it take forever)
        tar, tar_spec = saphires.io.read_vars( file_in['wavelength'].data, flux_cont_norm, header_df['file_token'].values[i_file], w_file = temp_ls_file_name, combine_all = False )

        # Make the saphires template spetrum structure. Read in the csv file from the code data directory and make into saphires structure
        template_csv  = pd.read_csv( os.path.join( config['paths']['code_dir'], 'data', config['radial_velocity']['template_file_name'] ) )
        template_spec = saphires.io.read_vars( template_csv['wavelength'], template_csv['flux'], config['radial_velocity']['template_file_name'], temp = True )

        # Prepare the spectra to run the broadening function
        tar_spec = saphires.utils.prepare( tar, tar_spec, template_spec, cr_trim = -0.2, oversample = 1 )

        ### Run broadening function calculation!
        
        try: # Put in try/except in case there are any issues with the BF computation
            tar_spec = saphires.bf.compute( tar, tar_spec, vel_width = config['radial_velocity']['bf_velocity_span'] )
        except RuntimeError:
            print( 'Issue with computing the broadening function!' )
            continue

        ### Now run analysis to smooth the broadening function
        
        try:
            tar_spec = saphires.bf.analysis( tar, tar_spec, sb = 'sb1', R = config['radial_velocity']['bf_smooth_res'], fit_trim = 20, text_out = False, prof = 'g' )
        except RuntimeError:
            print( 'Issue with fitting/smoothing the broadening function!' )
            continue
        
        # Remove the temporary ls file name
        os.remove( temp_ls_file_name )
        
        ### Bootstrap sample BF weight combination, to get RV and error in RV
        
        bootstrap_rv_samples = bootstrap_sample_bf_rvs( tar, tar_spec, config['radial_velocity']['n_bootstrap_samples'] )
        
        # Apply barycentric velocity correction to the samples
        bootstrap_rv_samples = bootstrap_rv_samples + barycorr_vel + np.median( bootstrap_rv_samples ) * barycorr_vel / constants.c.to('km/s').value

        # Check if all of the RV bootstrap samples are nans -- issue!
        if np.all( np.isnan( bootstrap_rv_samples ) ):
            print( 'ISSUE: All bootstrap RV samples are nans!' )
            continue

        # Plot the result of the RV sampling/BF combination
        plot_file_name = os.path.join( config['paths']['reduction_dir'], 'radial_velocity', 'bf_rv_result_{}.pdf'.format( header_df['file_token'].values[i_file] ) )
        plot_bootstrap_rv_result( tar, tar_spec, barycorr_vel, bootstrap_rv_samples, plot_file_name )
        
        # Calculate RV value and error
        bootstrap_rv_value = np.nanmedian( bootstrap_rv_samples )
        bootstrap_rv_error = median_abs_deviation( bootstrap_rv_samples, scale = 'normal', nan_policy = 'omit' )
    
        ### Append to output file!
        
        # Set up output BF arrays
        output_bf_arrays = np.full( ( 122, tar_spec[tar[0]]['vel'].size ), np.nan )
                
        for i_order, order in enumerate( orders_to_use['order'].values ):
            output_bf_arrays[order] = tar_spec[tar[i_order]]['bf_smooth']
        
        # Build a new HDU list with the radial velocity extension added. Appends to first 5 extensions in the input file (through continuum, just in case there is overwriting weirdness to avoid duplication)
        output_file = fits.HDUList( file_in[:5] + [ fits.ImageHDU( output_bf_arrays, name = 'radial velocity' ) ] )
        
        # Add the barycentric velocity correction to the header
        output_file[0].header['BARYCORR'] = ( barycorr_vel, 'Barycentric velocity correction (km/s)' )
        
        # Add radial velocity and error to the header, both primary header and the radial velocity extension
        for ext in [ 0, 5 ]:
            output_file[ext].header['RVBF']  = ( bootstrap_rv_value, 'Broadening function RV (km/s)' )
            output_file[ext].header['ERVBF'] = ( bootstrap_rv_error, 'Broadening function RV Error (km/s)' )
        
        # Add information for reconstructing the BF velocity array in the radial velocity extension header
        output_file['radial velocity'].header['VELSTART'] = ( tar_spec[tar[0]]['vel'][0], 'BF velocity array start (km/s)' )
        output_file['radial velocity'].header['VELSTEP']  = ( tar_spec[tar[0]]['vel_spacing'], 'BF velocity array spacing (km/s)' )
        output_file['radial velocity'].header['NVELPTS']  = ( tar_spec[tar[0]]['vel'].size, 'BF velocity array size' )
        
        # Add history to the primary header
        output_file[0].header['HISTORY'] = 'Radial velocity measured on {}'.format( datetime.strftime( datetime.now(), '%Y/%m/%d' ) )
        
        # Write out the file
        output_file.writeto( file_name, overwrite = True )

    return None


