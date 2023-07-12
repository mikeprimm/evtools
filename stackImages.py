import argparse
from genericpath import isfile
import os
import pathlib
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.time import Time
import numpy as np
import logging

import warnings
# UTC to BJD converter import
from pandas import isna
from reproject import reproject_interp

try:
    from .libs.stacks import buildMedianStack, scaleAndDemosaicImage
except ImportError:  # package import
    from libs.stacks import buildMedianStack, scaleAndDemosaicImage

warnings.simplefilter('ignore', category=FITSFixedWarning)

# create logger
logger = logging.getLogger('processExoplanetTransit')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

# Initialize parser
parser = argparse.ArgumentParser()
# Add input argument
parser.add_argument("-i", "--input", help = "Source image directory", required=True);
parser.add_argument("-d", "--darks", help = "Darks directory", required=True);
# Adding output argument
parser.add_argument("-o", "--output", help = "Output stacked image directory") 
parser.add_argument("--stacktime", help = "Number of seconds to stack (default 120)") 

# Read arguments from command line
try:
    args = parser.parse_args()
except argparse.ArgumentError:
    os.exit(1)
darksrcdir='darks'
if args.darks:
    darksrcdir = args.darks 
outputdir='output'
if args.output: 
    outputdir = args.output
inputdir='input'
if args.input:
    inputdir = args.input
stacktime = 120
if args.stacktime:
    stacktime = int(args.stacktime)
    
# Make output directory, if needed
pathlib.Path(outputdir).mkdir(parents=True, exist_ok=True)
tmppath = os.path.join(outputdir, "tmp")
pathlib.Path(tmppath).mkdir(parents=True, exist_ok=True)

ch2 = logging.FileHandler(os.path.join(outputdir, 'processExoplanetTransit.log'), encoding='utf-8', mode='w')
ch2.setLevel(logging.DEBUG)
ch2.setFormatter(formatter)
logger.addHandler(ch2)

darkfiles = []
# Go through the darks
for path in os.listdir(darksrcdir):
    dfile = os.path.join(darksrcdir, path)
    if (path.startswith('.')): continue
    # check if current path is a file
    if os.path.isfile(dfile):
        darkfiles.append(path)
darkfiles.sort()
# Build dark frame, if we have any to work with
dark = fits.HDUList()
if len(darkfiles) > 0:
    logger.info(f"Processing {len(darkfiles)} darks")
    dark = buildMedianStack(darksrcdir, darkfiles, "master-dark.fits")

# Go through the inputs
lightfiles = []
for path in os.listdir(inputdir):
    if (path.startswith('.')): continue
    dfile = os.path.join(inputdir, path)
    # check if current path is a file
    if os.path.isfile(dfile):
        lightfiles.append(path)
lightfiles.sort()

cnt = 0
stackedcnt = 0
timeaccumlist = []
timeaccumstart = 0
timeaccumra = 0
timeaccumdec = 0
mjdobs = 0
mjdend = 0

for f in lightfiles:
    try:
        lfile = os.path.join(inputdir, f)
        # Load file into list of HDU list 
        with fits.open(lfile) as hduList:
            imagedata = hduList[0].data
            print(f"imagedata.shape={imagedata.shape}")
            # Dark subtract, if we have darks
            if len(dark) > 0:
                # Clamp the data with the dark from below, so we can subtract without rollover
                dat = np.maximum(imagedata, dark[0].data)
                # And subtract the dark
                imagedata = np.subtract(dat, dark[0].data)
            print(f"imagedata.shape={imagedata.shape}")
            red, green, blue = scaleAndDemosaicImage(imagedata)
            tobs = hduList[0].header['MJD-OBS'] * 24 * 60 # MJD in minutes
            # End of accumulator?
            if (len(timeaccumlist) > 0) and ((timeaccumstart + stacktime) < tobs):
                stackout = os.path.join(outputdir, "stack-%04d.fits" % stackedcnt)
                tmpout = os.path.join(tmppath, "tmpstack.fits")
                timeaccumlist = []
                timeaccumstart = 0
            timeaccumlist.append(lfile)
            mjdend = hduList[0].header['MJD-END']
        cnt = cnt + 1
    except OSError as e:
        print("Error: file %s - %s (%s)" % (f, e.__class__, e))     
# Final accumulator?
if len(timeaccumlist) > 0:
    stackout = os.path.join(outputdir, "stack-%04d.fits" % stackedcnt)

print("Processed %d images into %d stacked images" % (cnt, stackedcnt))

