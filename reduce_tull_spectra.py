#!/usr/bin/env python

##### Imports #####

import argparse
import os
import yaml

import numpy as np

from astropy.io import fits

# Import pipeline and reduction modules

import tull_coude_reduction

from tull_coude_reduction.modules import ( 
    ccd_calibrations,
    image_processing,
    trace_echelle,
    extract_spectrum,
    wavelength_solve_and_calibrate,
    continuum_fit,
    radial_velocity,
    reduction_utils )

##### Set up for running the reduction pipeline! #####

### Read in the command line argument for the config file that defines the reduction run

parser = argparse.ArgumentParser( description = 'Tull coude reduction script' )
parser.add_argument( 'config_file', type = str, help = 'Configuration file with settings for the Tull coude reduction steps' )
args = parser.parse_args()

# Read in the YAML config file for this reduction run
config_file = yaml.safe_load( open( args.config_file, 'r' ) )

print( 'PIPELINE START: Running reduction pipeline for night at data path: {}'.format( config_file['paths']['working_dir'] ) )

### Add an entry to the config paths section pointing to the data folder in the code directory

config_file['paths']['code_data_dir'] = os.path.join( tull_coude_reduction.__path__[0], 'data' )

### Set up the reduction directories for holding the pipeline output

# cd into the working directory
os.chdir( config_file['paths']['working_dir'] )

# Make the directory tree required within the working directory
for dir_name in config_file['paths']['sub_dir_list']:
    os.makedirs( os.path.join( config_file['paths']['reduction_dir'], dir_name ), exist_ok = True )

### Read in the FITS files and pull header information

# Write the header information to a CSV file
header_info = reduction_utils.make_header_manifest( config_file['general']['header_info_file_name'] )

# Also generate a list of the object names from the file headers in lower case, for ease of searching through them later
object_names_lowercase = np.array( [ object_name.lower() for object_name in header_info['object'].values ] )

##### Now run the reduction pipeline modules! #####

### Check the config general "do_all_steps" flag -- if True, override individual module "do_steps" and do everything!

if config_file['general']['do_all_steps']:
    
    # Go through each of the config sections
    for config_key in config_file.keys():
        
        # Only work on sections that have a "do_step" key
        if 'do_step' in config_file[config_key]:
            config_file[config_key]['do_step'] = True
        

### Build CCD calibration files -- bias, flat, and bad pixel mask

# Get the file indices for the bias frames in the header info file, with the image type string defined in the config file (Tull coude: 'zero')
bias_frame_indices = np.where( header_info['image_type'].values == config_file['calibrations']['bias_image_type'] )[0]

# Get the file indices fro the flat frames in the header info file
# The config file defines: 1) the image type string for flats (Tull coude: 'flat') and 2) any object names to discard
flat_frame_indices = np.where( ( header_info['image_type'].values == config_file['calibrations']['flat_image_type'] ) & 
                               ( [ object_name not in config_file['calibrations']['flat_object_names_to_discard'] for object_name in header_info['object'].values ] ) )[0]

# Run the CCD calibration file making step if config file says to do so
if config_file['calibrations']['do_step']:
    print( 'MODULE START: Making CCD calibration files.' )

    # Run CCD calibration file making module
    ccd_calibrations.build_calibrations( header_info, bias_frame_indices, flat_frame_indices, config_file )
    
### Now process and generate the 2D echellogram images for the files to be reduced, using the CCD calibration files

# Process and produce images if config file says to do so
if config_file['image_process']['do_step']:
    print( 'MODULE START: Image processing.' )

    # Get the file indices for the arc lamp frames to process, to exclude any testing/focusing frames
    arc_frame_indices = np.where( np.logical_and( header_info['image_type'].values == config_file['image_process']['arc_lamp_image_type'], np.any( [ object_names_lowercase == name for name in config_file['image_process']['valid_arc_lamp_object_names'] ], axis = 0 ) ) )[0]
    
    # Get the file indices for the science target frames to process, ecluding any test frames
    object_frame_indices = np.where( np.logical_and( header_info['image_type'].values == 'object', np.all( [ object_names_lowercase != name for name in config_file['image_process']['invalid_science_object_names'] ], axis = 0 ) ) )[0]

    # Concatenate the two together!
    frame_indices_to_process = np.sort( np.concatenate( [ arc_frame_indices, object_frame_indices ] ) )
    
    # Read in the CCD calibration files necessary for image processing
    super_bias     = fits.open( os.path.join( config_file['paths']['reduction_dir'], 'cals', 'super_bias.fits' ) )
    flat_field     = fits.open( os.path.join( config_file['paths']['reduction_dir'], 'cals', 'flat_field.fits' ) )
    bad_pixel_mask = fits.open( os.path.join( config_file['paths']['reduction_dir'], 'cals', 'bad_pixel_mask.fits' ) )
    
    # Run image processing module
    image_processing.build_images( frame_indices_to_process, super_bias, flat_field, bad_pixel_mask, header_info, config_file )
    
### Now get the trace for the echelle orders

# Find and fit the echellogram traces if the config file says to
if config_file['trace']['do_step']:
    print( 'MODULE START: Tracing echelle orders.' )
    
    # Read in the flat field
    flat_field = fits.open( os.path.join( config_file['paths']['reduction_dir'], 'cals', 'flat_field.fits' ) )
    
    # Run trace finding module
    trace_echelle.get_trace( flat_field, config_file )

### Extract 1D spectra from the processed 2D images

# If the config file says to
if config_file['extraction']['do_step']:
    print( 'MODULE START: Extracting 1D spectra.' )
    
    # Read in the trace file
    trace_file = fits.open( os.path.join( config_file['paths']['reduction_dir'], 'trace', 'trace.fits' ) )
        
    ## Extract the arc lamp spectra -- separate from the science object spectra because no optimal extraction
    
    # The file indices for the arc lamp frames to extract
    arc_frame_indices = np.where( np.logical_and( header_info['image_type'].values == config_file['image_process']['arc_lamp_image_type'], np.any( [ object_names_lowercase == name for name in config_file['image_process']['valid_arc_lamp_object_names'] ], axis = 0 ) ) )[0]

    # Run the extraction for arc lamp frames
    extract_spectrum.extract_spectrum( arc_frame_indices, trace_file[2].data, header_info, config_file['extraction']['lamp_extract_type'], config_file['extraction']['lamp_background_subtract'], config_file )
    
    ## Extract the science spectra
    
    # Leave this here for now -- want to add regex for excluding solar port data
    not_object_names = [ 'test', 'solar', 'sol port', 'solar port', 'solport', 'solarport', 'solar port halpha' ]

    object_frame_indices = np.where( np.logical_and( header_info['image_type'].values == 'object', np.all( [ object_names_lowercase != name for name in not_object_names ], axis = 0 ) ) )[0]

    # Run the extraction module for the science frames    
    extract_spectrum.extract_spectrum( object_frame_indices, trace_file[2].data, header_info, config_file['extraction']['science_extract_type'], config_file['extraction']['science_background_subtract'], config_file )
    
### Make the wavelength solution and calibrate the spectra

if config_file['wavecal']['do_step']:
    print( 'MODULE START: Wavelength calibration.' )
    
    ## Get the file indices of the arc lamp frames to use for wavelength solution -- which allows for a minimum exposure time set in the config file
    
    # First the indices where the type and object names are for arc lamps
    arc_frame_indices = np.where( np.logical_and( header_info['image_type'].values == config_file['image_process']['arc_lamp_image_type'], np.any( [ object_names_lowercase == name for name in config_file['image_process']['valid_arc_lamp_object_names'] ], axis = 0 ) ) )[0]
    
    # And now apply the minimum exposure time criterion
    arc_frame_indices = np.intersect1d( arc_frame_indices, np.where( header_info['exp_time'].values >= config_file['wavecal']['min_arc_exp_time'] )[0] )
    
    # Run the wavelength solution module
    wavelength_solve_and_calibrate.wavelength_solution( arc_frame_indices, header_info, config_file )

    ## Get the file indices for the science spectra to wavelength calibrate and calibrate them
    
    # Just take all of the files that have extracted spectra files other than the arc lamp frames used for the wavelength solution
    frame_indices_to_wavecal = [ i for i in range( header_info.shape[0] ) if os.path.exists( os.path.join( config_file['paths']['reduction_dir'], 'spectrum_files/tullcoude_{}_spectrum.fits'.format( header_info['file_token'].values[i] ) ) ) ]
    frame_indices_to_wavecal = np.setdiff1d( frame_indices_to_wavecal, arc_frame_indices )
    
    # Run the wavelength calibration module
    wavelength_solve_and_calibrate.wavelength_calibrate( frame_indices_to_wavecal, arc_frame_indices, header_info, config_file )
    
### Continuum fit the extracted science spectra

if config_file['continuum_fit']['do_step']:
    print( 'MODULE START: Fitting continuum' )
    
    # The file indices of the science frames to continuum fit
    not_object_names = [ 'test', 'solar', 'sol port', 'solar port', 'solport', 'solarport', 'solar port halpha' ]

    object_frame_indices = np.where( np.logical_and( header_info['image_type'].values == 'object', np.all( [ object_names_lowercase != name for name in not_object_names ], axis = 0 ) ) )[0]
    
    # Run the spline continuum fitting module
    continuum_fit.fit_spectra_continuum( object_frame_indices, header_info, config_file )
    
if config_file['radial_velocity']['do_step']:
    print( 'MODULE START: Calculating radial velocities' )
    
    # The file indices of the science frames to measure radial velocities for
    not_object_names = [ 'test', 'solar', 'sol port', 'solar port', 'solport', 'solarport', 'solar port halpha' ]

    object_frame_indices = np.where( np.logical_and( header_info['image_type'].values == 'object', np.all( [ object_names_lowercase != name for name in not_object_names ], axis = 0 ) ) )[0]

    # Run the BF computation and RV measurement module
    radial_velocity.measure_radial_velocity( object_frame_indices, header_info, config_file )
    
    # Run function to generate the night's compiled RV information csv
    csv_file_name = os.path.join( config_file['paths']['reduction_dir'], 'radial_velocity', 'compiled_rv_info.csv' )
    radial_velocity.make_rv_compiled_excel( object_frame_indices, csv_file_name, header_info, config_file )
        
print( 'Everything is done.' )
    
