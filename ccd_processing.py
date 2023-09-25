""" Library of functions for the CCD processing steps to reduce data from the Tull coude spectrograph
Last updated: DMK, 9/24/2023
"""

##### Imports

import numpy as np

from astropy.io import fits

##### Functions

def build_super_bias( bias_file_names, read_noise ):
    """ Function to read in the bias images and combine them into a super bias to apply to other frames
    
    Parameters
    ----------
    bias_file_names : list
        List of file names (type str) for the raw bias frames
        
    read_noise : float
        The read noise (in X units)

    Returns
    -------
    super_bias : HDUList
        Multi-extension FITS ImageHDU list with the median combined super bias and associated error

    """
    
    # Initialize an array to hold the bias images
    bias_cube = []
    
    for i_file, file_name in enumerate( bias_file_names ):
        
        # Read in the bias file
        bias_image = fits.getdata( file_name )
        
        # Append to the cube of all bias images
        bias_cube.append( bias_image )
        
    # Get the super bias -- median combine the bias images, and add in read noise to bias error
    super_bias_flux  = np.nanmedian( bias_cube, axis = 0 )
    super_bias_error = np.sqrt( super_bias_flux + read_noise ** 2.0 )
    
    ### Output -- into FITS object
    
    # Put the super bias data into FITS ImageHDU objects
    super_bias_flux_hdu  = fits.ImageHDU( super_bias_flux, name = 'bias flux' )
    super_bias_error_hdu = fits.ImageHDU( super_bias_error, name = 'bias error' )
    
    # Output to a FITS multi-extension HDU object (with an empty Primary HDU to hold some header information)
    super_bias = fits.HDUList( [ fits.PrimaryHDU(), super_bias_flux_hdu, super_bias_error_hdu ] )
    
    # Add some keywords to the output super bias object
    super_bias[0].header['nbiases'] = ( len( bias_file_names ), 'Number of bias frames combined' )
    super_bias[0].header['rdnoise'] = ( read_noise, 'Read noise added to error' )
    super_bias[0].header['history'] = 'Median combined super bias frame'
        
    return super_bias

def build_flat_field( flat_file_names, read_noise, super_bias ):
    """ Function to read in flat field images and combine them into the flat field for image processing

    Parameters
    ----------
    flat_file_names : list
        List of file names (type str) for the raw flat field frames
        
    read_noise : float
        The read noise (in X units)
        
    super_bias : HDUList
        Multi-extension FITS ImageHDU list with the median combined super bias and associated error created with build_super_bias.

    Returns
    -------
    flat_field : HDUList
        Multi-extension FITS ImageHDU list with the median combined flat field and associated error.

    """
    
    # Initialize an array to hold the flat field images
    flat_cube = []
    
    for i_file, file_name in enumerate( flat_file_names ):
        
        # Read in the bias file
        flat_image = fits.getdata( file_name )
        
        # Append to the cube of all bias images
        flat_cube.append( flat_image )
    
    # Get the flat field -- median combine the flat field images and subtract the super bias
    flat_field_flux  = np.nanmedian( flat_cube, axis = 0 )
    flat_field_flux -= super_bias['bias flux'].data
    
    # Get the flat field errors (flat field, read noise, and super bias)
    flat_field_error = np.sqrt( flat_field_flux + read_noise ** 2.0 + super_bias['bias error'].data ** 2.0 )
    
    ### Output
    
    # Put the flat field flux and error into ImageHDU objects, with minor math to make the flat field scale between 0 and 1
    flat_field_flux_hdu  = fits.ImageHDU( ( flat_field_flux - np.nanmin( flat_field_flux ) ) / np.nanmax( flat_field_flux ), name = 'flat flux' )
    flat_field_error_hdu = fits.ImageHDU( flat_field_error / np.nanmax( flat_field_flux - np.nanmin( flat_field_flux ) ), name = 'flat error' )
    
    # Put into an HDUList
    flat_field = fits.HDUList( [ fits.PrimaryHDU(), flat_field_flux_hdu, flat_field_error_hdu ] )
        
    # Add some keywords to the output super bias object
    flat_field[0].header['nflats']  = ( len( flat_file_names ), 'Number of flat field frames combined' )
    flat_field[0].header['rdnoise'] = ( read_noise, 'Read noise added to error' )
    flat_field[0].header['history'] = 'Median combined flat field frame'
    flat_field[0].header['history'] = 'Super bias subtracted'

    return flat_field

def make_bad_pixel_mask( super_bias, flat_field, bias_bpm_percentile, flat_field_bpm_limit ):
    """

    Parameters
    ----------
    super_bias : TYPE
        DESCRIPTION.
        
    flat_field : TYPE
        DESCRIPTION.
        
    bias_bpm_percentile : TYPE
        DESCRIPTION.

    flat_field_bpm_limit : TYPE
        DESCRIPTION.

    Returns
    -------
    bad_pixel_mask : TYPE
        DESCRIPTION.

    """
    
    # Get the super bias flux value for the percentile above which to label bad pixels, as defined in the top level config
    bias_percentile_cut_value = np.nanpercentile( super_bias['bias flux'].data, bias_bpm_percentile )
    
    # Define the bad pixel mask as indices where the super bias flux is above the percentile limit, or the flat field is very small
    bad_pixel_mask = np.where( ( super_bias['bias flux'].data > bias_percentile_cut_value ) | ( flat_field['flat flux'].data < flat_field_bpm_limit ) )
    
    return bad_pixel_mask

##### Wrapper for building all calibration files

def build_calibrations( header_df, bias_frame_indices, flat_frame_indices, config ):
    """
    

    Parameters
    ----------
    header_df : TYPE
        DESCRIPTION.
    bias_frame_indices : TYPE
        DESCRIPTION.
    flat_frame_indices : TYPE
        DESCRIPTION.
    config : TYPE
        DESCRIPTION.

    Returns
    -------
    super_bias : TYPE
        DESCRIPTION.
    flat_field : TYPE
        DESCRIPTION.
    bad_pixel_mask : TYPE
        DESCRIPTION.

    """
    return super_bias, flat_field, bad_pixel_mask