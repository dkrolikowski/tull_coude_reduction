""" Functions to measure radial velocities from extracted echelle spectra.

Created by DMK on 11/20/2023
Last updated by DMK on 11/20/2023
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
import saphires

import pdb

from radec_conv import MS2Deg

# import numpy as np
# import pandas as pd

# import saphires as saph

# import glob
# import pickle
# import tqdm


##### Functions

def make_saphires_ls_file( orders_to_use, data_wavelength ):
    
    output_ls_file = open( 'temp.ls', 'w' )
    
    for i_order, order in enumerate( orders_to_use['order'].values ):
        
        wave_start, wave_end = data_wavelength[order].min() - 10, data_wavelength[order].max() + 10
        
        if orders_to_use['flag'].values[i_order] == 'start':
            wave_start = orders_to_use['wave_start'].values[i_order]
        elif orders_to_use['flag'].values[i_order] == 'end':
            wave_end = orders_to_use['wave_end'].values[i_order]
        
        output_ls_file.write( '{} {:.2f}-{:.2f}\n'.format( order, wave_start, wave_end ) )
        
    output_ls_file.close()
        
    return None

def bootstrap_sample_bf_rvs( bf_tar, bf_tar_spec, num_samples ):

    # Do bootstrap BF combination for getting the RV value!
    rv_samples = np.zeros( num_samples )

    for i_sample in range( num_samples ):

        orders_to_sample = np.random.randint( 0, len( bf_tar ), len( bf_tar ) )

        v_comb, bf_comb, _, _ = saphires.bf.weight_combine( bf_tar[orders_to_sample], bf_tar_spec, std_perc = 0.1 )

        p0 = [ np.max( bf_comb ), v_comb[np.argmax( bf_comb )], 1, np.median( bf_comb ) ]

        # Now fit
        try:
            fit, errs = curve_fit( saphires.utils.gaussian_off, v_comb, bf_comb, p0 = p0 )
        except RuntimeError:
            continue

        rv_samples[i_sample] = fit[1] + bf_tar_spec[bf_tar[0]]['rv_shift']

    return rv_samples

def plot_bootstrap_rv_result( bf_tar, bf_tar_spec, bc_vel, bootstrap_rv_samples, bootstrap_rv_value, bootstrap_rv_error, file_name ):
    
    # Combine the BFs
    velocity_arr, combined_bf, _, _ = saphires.bf.weight_combine( bf_tar, bf_tar_spec, std_perc = 0.1 )
    
    # Make the figure -- two panels
    fig, axs = plt.subplots( 1, 2, num = 1, clear = True )
    fig.set_size_inches( 15, 7 )
    
    ### Left panel: BF and Gaussian fit
    
    # Plot the combined BF
    axs[0].plot( velocity_arr + bc_vel, combined_bf, c = '#323232', lw = 1.5, label = None, zorder = 3 )
    
    # Plot the Gaussian fit!
    fit, errs = curve_fit( saphires.utils.gaussian_off, velocity_arr, combined_bf, p0 = [ np.max( combined_bf ), velocity_arr[np.argmax( combined_bf )], 1, np.median( combined_bf ) ] )
    fit_rv    = fit[1] + bc_vel + fit[1] * bc_vel / constants.c.to('km/s').value
    
    axs[0].plot( velocity_arr + bc_vel, saphires.utils.gaussian_off( velocity_arr, *fit ), '--', c = '#bf3465', lw = 1, label = 'Fit RV: {:.3f}'.format( fit_rv ), zorder = 3 )
    
    # Plot vertical lines for the bootstrap sample RVs
    axs[0].axvline( x = bootstrap_rv_value, c = '#dfa5e5', label = 'Bootstrap RV: {:.3f} +/- {:.3f}'.format( bootstrap_rv_value, bootstrap_rv_error ), lw = 1.5, zorder = 2 )
    axs[0].axvspan( bootstrap_rv_value - 3 * bootstrap_rv_error, bootstrap_rv_value + 3 * bootstrap_rv_error, color = '#dfa5e5', alpha = 0.1, zorder = 1 )

    axs[0].set_xlabel( 'Velocity (km/s)' )
    axs[0].set_ylabel( 'Brodening Function' )
    axs[0].legend( fontsize = 'small' )
    
    ### Right panel: RV bootstrap samples
    
    axs[1].hist( bootstrap_rv_samples, bins = 50, histtype = 'step', color = '#323232', lw = 1.5, range = ( bootstrap_rv_value - 10 * bootstrap_rv_error, bootstrap_rv_value + 10 * bootstrap_rv_error ), zorder = 3 )
    
    # Plot vertical lines for the RV values
    axs[1].axvline( x = bootstrap_rv_value, c = '#dfa5e5', label = 'Bootstrap RV: {:.3f} +/- {:.3f}'.format( bootstrap_rv_value, bootstrap_rv_error ), lw = 1.5, zorder = 2 )
    axs[1].axvspan( bootstrap_rv_value - bootstrap_rv_error, bootstrap_rv_value + bootstrap_rv_error, color = '#dfa5e5', alpha = 0.1, zorder = 1 )

    axs[1].set_xlabel( 'Bootstrap Sample RV (km/s)' )
    axs[1].set_ylabel( 'Histogram Counts' )
    axs[1].legend( fontsize = 'small' )

    plt.savefig( file_name, bbox_inches = 'tight', pad_inches = 0.05 )
    
    return None

##### Main script

def measure_radial_velocity( file_indices, header_df, config ):
    
    ### Loop through each of the files to measure RV for
    for i_file in file_indices:
        
        ### Data preparation -- read in, normalize, get barycentric velocity correction
        
        # Read in the spectra file
        file_name = os.path.join( config['paths']['reduction_dir'], 'spectrum_files', 'tullcoude_{}_spectrum.fits'.format( header_df['file_token'].values[i_file] ) )
        file_in   = fits.open( file_name )
                
        # Continuum normalize the spectrum
        flux_cont_norm = file_in['extracted flux'].data / file_in['continuum'].data

        # Get the data midtime. If flux weighted midtime is in the header use that, otherwise the header OBSJD + half exposure time
        if 'emfwmtim' in file_in[0].header:
            obs_mid_jd = Time( file_in[0].header['DATE-OBS'] + 'T' + file_in[0].header['emfwmtim'], format = 'isot' ).jd
        else:
            obs_mid_jd = header_df['obs_jd'].values[i_file] + 0.5 * header_df['exp_time'].values[i_file] / 60. / 60. / 24.
                
        ## BC velocity
        
        # Get RA/Dec from header into degrees
        sky_coord = SkyCoord( file_in[0].header['ra'], file_in[0].header['dec'], unit = ( u.hourangle, u.deg ) )
                
        star_info    = { 'ra': sky_coord.ra.value, 'dec': sky_coord.dec.value, 'epoch': 2451545.0 }
        barycorr_vel = barycorrpy.get_BC_vel( obs_mid_jd, obsname = 'McDonald', leap_update = False, **star_info )[0][0] / 1000.0
        
        ### Saphires broadening function prep -- create ls file, make data and template structures
        
        # Read in the pre-defined orders to use file
        orders_to_use = pd.read_csv( os.path.join( config['paths']['code_dir'], 'data', config['radial_velocity']['orders_to_use_file_name'] ) )
        
        # Make the ls file
        make_saphires_ls_file( orders_to_use, file_in['wavelength'].data )
        
        # Make the saphires target data structures. Don't combine all orders (unnecessary and makes it take forever)
        tar, tar_spec = saphires.io.read_vars( file_in['wavelength'].data, flux_cont_norm, header_df['file_token'].values[i_file], w_file = 'temp.ls', combine_all = False )

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
            print( 'Issue with fitting the broadening function!' )
            continue
        
        ### Bootstrap sample BF weight combination, to get RV and error in RV
        
        bootstrap_rv_samples = bootstrap_sample_bf_rvs( tar, tar_spec, config['radial_velocity']['n_bootstrap_samples'] )
        
        # Apply barycentric velocity correction to the samples
        bootstrap_rv_samples = bootstrap_rv_samples + barycorr_vel + np.median( bootstrap_rv_samples ) * barycorr_vel / constants.c.to('km/s').value
        
        # Calculate RV value and error
        rv_value = np.median( bootstrap_rv_samples )
        rv_error = median_abs_deviation( bootstrap_rv_samples, scale = 'normal' )
        
        plot_file_name = os.path.join( config['paths']['reduction_dir'], 'radial_velocity', 'bf_rv_result_{}.pdf'.format( header_df['file_token'].values[i_file] ) )
        plot_bootstrap_rv_result( tar, tar_spec, barycorr_vel, bootstrap_rv_samples, rv_value, rv_error, plot_file_name )
        
        ### Append to output file!
        
        # Set up output BF arrays
        output_bf_arrays = np.full( ( 122, tar_spec[tar[0]]['vel'].size ), np.nan )
                
        for i_order, order in enumerate( orders_to_use['order'].values ):
            output_bf_arrays[order] = tar_spec[tar[i_order]]['bf_smooth']
        
        # Build a new HDU list with the continuum extension added
        output_file = fits.HDUList( file_in[:5] + [ fits.ImageHDU( output_bf_arrays, name = 'radial velocity' ) ] )
        
        # Add radial velocity and error to the header
        output_file[0].header['RVBF']  = ( rv_value, 'Broadening function RV (km/s)' )
        output_file[0].header['ERVBF'] = ( rv_error, 'Broadening function RV Error (km/s)' )
        
        # Add information for reconstructing the BF velocity array in the radial velocity extension header
        output_file['radial velocity'].header['VELSTART'] = ( tar_spec[tar[0]]['vel'][0], 'BF velocity array start (km/s)' )
        output_file['radial velocity'].header['VELSTEP']  = ( tar_spec[tar[0]]['vel_spacing'], 'BF velocity array spacing (km/s)' )
        output_file['radial velocity'].header['NVELPTS']  = ( tar_spec[tar[0]]['vel'].size, 'BF velocity array size' )
        
        # Add history to the primary header
        output_file[0].header['HISTORY'] = 'Radial velocity measured on {}'.format( datetime.strftime( datetime.now(), '%Y/%m/%d' ) )
        
        # Write out the file
        output_file.writeto( file_name, overwrite = True )

    return None






