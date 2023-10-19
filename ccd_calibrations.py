""" Library of functions to build the CCD calibration files for image processing to reduce data from the Tull coude spectrograph

Created by DMK on 9/24/2023
Last updated by DMK on 10/14/2023
"""

##### Imports

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from os import path

##### Functions

def build_super_bias( bias_file_names, read_noise ):
    """ Function to read in the bias images and median combine them into a super bias to apply to other frames
    
    Parameters
    ----------
    bias_file_names : list
        List of file names (type str) for the raw bias frames
        
    read_noise : float
        The read noise (in ADU, converted from the header value which is given in electrons)

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
    """ Function to read in flat field images and median combine them into the flat field for image processing

    Parameters
    ----------
    flat_file_names : list
        List of file names (type str) for the raw flat field frames
        
    read_noise : float
        The read noise (in ADU, converted from the header value which is given in electrons)
        
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
    flat_field_median = np.nanmedian( flat_cube, axis = 0 )
    flat_field_flux   = flat_field_median - super_bias['bias flux'].data
    
    # Get the flat field errors (flat field, read noise, and super bias)
    flat_field_error = np.sqrt( flat_field_median + read_noise ** 2.0 + super_bias['bias error'].data ** 2.0 )
    
    ### Output
    
    # Put the flat field flux and error into ImageHDU objects, with minor math to make the flat field scale between 0 and 1
    flat_field_flux_hdu  = fits.ImageHDU( ( flat_field_flux - np.nanmin( flat_field_flux ) ) / np.nanmax( flat_field_flux - np.nanmin( flat_field_flux ) ), name = 'flat flux' )
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
    """ Function to create a bad pixel mask using the super bias and flat field.

    Parameters
    ----------
    super_bias : HDUList
        Multi-extension FITS ImageHDU list with the median combined super bias and associated error created with build_super_bias.
        
    flat_field : HDUList
        Multi-extension FITS ImageHDU list with the median combined flat field and associated error created with build_flat_field.
        
    bias_bpm_percentile : float
        The percentile above which super bias counts are marked as bad pixels (hot pixels). This is defined in the config file.

    flat_field_bpm_limit : float
        The flat field value below which a pixel is marked as bad. This is defined in the config file.

    Returns
    -------
    bad_pixel_mask : FITSHDU
        The bad pixel mask created using the super bias and flat field. Bad pixels are valued 0, good pixels are valued 1.

    """
    
    # Get the super bias flux value for the percentile above which to label bad pixels, as defined in the top level config
    bias_percentile_cut_value = np.nanpercentile( super_bias['bias flux'].data, bias_bpm_percentile )
    
    # Define the bad pixel mask as indices where the super bias flux is above the percentile limit, or the flat field is very small
    bad_pixel_mask_indices = np.where( ( super_bias['bias flux'].data > bias_percentile_cut_value ) | ( flat_field['flat flux'].data < flat_field_bpm_limit ) )
    
    ### Output the bad pixel mask as a fits file: 0s are bad pixels, 1s are good pixels
    
    bad_pixel_mask = np.ones( super_bias['bias flux'].data.shape )
    
    bad_pixel_mask[bad_pixel_mask_indices] = 0
    
    # FITS HDU
    bad_pixel_mask_hdu = fits.PrimaryHDU( bad_pixel_mask )
    
    bad_pixel_mask_hdu.header['history'] = 'Bias above {:.2f}'.format( bias_percentile_cut_value )
    bad_pixel_mask_hdu.header['history'] = 'Flat field below {:.5f}'.format( flat_field_bpm_limit )
    
    return bad_pixel_mask_hdu

def cal_image_2d_plot( image, figure_title, file_name, bpm = None ):
    """ Function to make a 2D image plot of an input 2D data array. Generalized to plot any input image.

    Parameters
    ----------
    image : array
        The 2D array to show.
        
    title : str
        The title of the plot (here to mark what type of image is plotted).
        
    file_name : str
        The file name to save the plot to.
        
    bpm : array, optional
        The 2D array defining the bad pixel mask. This can be input and overplotted the input data image. The default is None.

    Returns
    -------
    None.

    """
    
    plt.clf()
    
    # Plot the image
    plt.imshow( np.log10( image ), cmap = 'gray', aspect = 'auto', interpolation = 'none' )
    
    # If the bad pixel mask is given, plot the bad pixels on top!
    if bpm is not None:
        bpm_pixels = np.where( bpm == 0 )
        
        # Invert the x and y BPM for imshow
        plt.plot( bpm_pixels[1], bpm_pixels[0], '.', ms = 1, c = '#db6d1b', label = 'Bad pixel mask' )
        
        plt.legend( fontsize = 'small' )
        
    # Add colorbar
    plt.colorbar()
    
    # Add title passed to the function
    plt.title( figure_title )
    
    # Save the figure
    plt.savefig( file_name, bbox_inches = 'tight', pad_inches = 0.05 )
    plt.clf()
        
    return None

##### Wrapper for building all calibration files

def build_calibrations( header_df, bias_frame_indices, flat_frame_indices, config ):
    """ Main function to build all of the CCD calibrations, using the above functions for the individual calibrations.

    Parameters
    ----------
    header_df : pandas DataFrame
        The compiled information from the file headers.
        
    bias_frame_indices : list
        List of indices of the bias frames for the header dataframe.
        
    flat_frame_indices : list
        List of indices of the flat frames for the header dataframe.
        
    config : dict
        The overall config file defined using YAML with all parameters for running the reduction and analysis pipeline.

    Returns
    -------
    super_bias : HDUList
        Multi-extension FITS ImageHDU list with the median combined super bias and associated error created with build_super_bias.
        
    flat_field : HDUList
        Multi-extension FITS ImageHDU list with the median combined flat field and associated error created with build_flat_field.
        
    bad_pixel_mask : FITSHDU
        The bad pixel mask created with make_bad_pixel_mask using the super bias and flat field.

    """
    
    ### Build the files
    
    # First the super bias
    print( 'Reading bias files and creating super bias' )
    
    super_bias = build_super_bias( header_df['file_name'].values[bias_frame_indices], ( header_df['read_noise'] / header_df['gain'] ).values[bias_frame_indices][0] )
    super_bias.writeto( path.join( config['paths']['reduction_dir'], 'cals', 'super_bias.fits' ), overwrite = True, checksum = True )
    
    # Next the flat field
    print( 'Reading flat files and creating flat field' )
    
    flat_field = build_flat_field( header_df['file_name'].values[flat_frame_indices], ( header_df['read_noise'] / header_df['gain'] ).values[flat_frame_indices][0], super_bias )
    flat_field.writeto( path.join( config['paths']['reduction_dir'], 'cals', 'flat_field.fits' ), overwrite = True, checksum = True )
    
    # Last the bad pixel mask
    print( 'Building the bad pixel mask from the super bias and flat field' )
    
    bad_pixel_mask = make_bad_pixel_mask( super_bias, flat_field, config['calibrations']['bias_bpm_percentile'], config['calibrations']['flat_field_bpm_limit'] )
    bad_pixel_mask.writeto( path.join( config['paths']['reduction_dir'], 'cals', 'bad_pixel_mask.fits' ), overwrite = True, checksum = True )
    
    ### Make quick check plots for the calibration files!
    
    # First for the super bias
    super_bias_title = 'Super Bias, Median Value = {0:.2f}'.format( np.nanmedian( super_bias['bias flux'].data ) )
    cal_image_2d_plot( super_bias['bias flux'].data, super_bias_title, path.join( config['paths']['reduction_dir'], 'cals', 'super_bias.pdf' ) )
    
    # Next for the flat fiel
    flat_field_title = 'Flat Field'
    cal_image_2d_plot( flat_field['flat flux'].data, flat_field_title, path.join( config['paths']['reduction_dir'], 'cals', 'flat_field.pdf' ) )
    
    # Last for the bad pixel mask
    bpm_title = 'Bad Pixel Mask, plotted on Super Bias'
    cal_image_2d_plot( super_bias['bias flux'].data, bpm_title, path.join( config['paths']['reduction_dir'], 'cals', 'bad_pixel_mask.pdf' ), bpm = bad_pixel_mask.data )
    
    return None

