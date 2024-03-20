import argparse
from genericpath import isfile
import os
import pathlib
import shutil
import logging
from tkinter import E
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from astropy.wcs.wcs import NoConvergence
import astropy.units as u
import numpy as np
import subprocess
import warnings
from colour_demosaicing import demosaicing_CFA_Bayer_bilinear

try:
    from .libs.stacks import buildMedianStack, buildMasterFlatStack
except ImportError:  # package import
    from libs.stacks import buildMedianStack, buildMasterFlatStack
try:
    from .libs.exofop import exofop_getcompositeinfo, exofop_getticid
except ImportError: 
    from libs.exofop import exofop_getcompositeinfo, exofop_getticid

# UTC to BJD converter import
from barycorrpy import utc_tdb
from pandas import isna

warnings.simplefilter('ignore', category=FITSFixedWarning)

# create logger
logger = logging.getLogger('calibrateImages')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

# Initialize parser
parser = argparse.ArgumentParser()
# Add input argument
parser.add_argument("-d", "--darks", help = "Dark files source directory");
parser.add_argument("-s", "--science", help = "Science files source directory");
parser.add_argument("-df", "--darkflats", help = "Dark flat files source directory");
parser.add_argument("-f", "--flats", help = "Flat files source directory");
# Adding output argument
parser.add_argument("-o", "--output", help = "Output directory") 

# Read arguments from command line
try:
    args = parser.parse_args()
except argparse.ArgumentError:
    os.exit(1)
outputdir='output'
if args.output: 
    outputdir = args.output
# Add file logger
pathlib.Path(outputdir).mkdir(parents=True, exist_ok=True)
ch2 = logging.FileHandler(os.path.join(outputdir, 'calibrateImages.log'), encoding='utf-8', mode='w')
ch2.setLevel(logging.DEBUG)
ch2.setFormatter(formatter)
logger.addHandler(ch2)

logger.info("args=%s" % args)

darksrcdir='darks'
if args.darks:
    darksrcdir = args.darks 
sciencesrcdir='science'
if args.science:
    sciencesrcdir = args.science
darkflatsrcdir=None
if args.darkflats:
    darkflatsrcdir = args.darkflats
flatsrcdir=None
if args.flats:
    flatsrcdir = args.flats

# Make output directory, if needed
pathlib.Path(outputdir).mkdir(parents=True, exist_ok=True)
darkpath = os.path.join(outputdir, "darks")
pathlib.Path(darkpath).mkdir(parents=True, exist_ok=True)
sciencepath = os.path.join(outputdir, "science")
pathlib.Path(sciencepath).mkdir(parents=True, exist_ok=True)
badsciencepath = os.path.join(outputdir, "science-rej")
pathlib.Path(badsciencepath).mkdir(parents=True, exist_ok=True)
tmppath = os.path.join(outputdir, "tmp")
pathlib.Path(tmppath).mkdir(parents=True, exist_ok=True)
if args.darkflats:
   darkflatpath = os.path.join(outputdir, "darkflats")
   pathlib.Path(darkflatpath).mkdir(parents=True, exist_ok=True)
if args.flats:
   flatpath = os.path.join(outputdir, "flats")
   pathlib.Path(flatpath).mkdir(parents=True, exist_ok=True)

darkfiles = []
calstat = ""
# Go through the darks
for path in os.listdir(darksrcdir):
    dfile = os.path.join(darksrcdir, path)
    if (path.startswith('.')): continue
    if (path.startswith('master-dark.fits')): continue
    # check if current path is a file
    if os.path.isfile(dfile):
        darkfiles.append(path)
darkfiles.sort()
# Go through the lights
lightfiles = []
for path in os.listdir(sciencesrcdir):
    if (path.startswith('.')): continue
    dfile = os.path.join(sciencesrcdir, path)
    # check if current path is a file
    if os.path.isfile(dfile):
        lightfiles.append(path)
lightfiles.sort()
# Do dark flats, if any
darkflatfiles = []
if darkflatsrcdir:
    for path in os.listdir(darkflatsrcdir):
        dfile = os.path.join(darkflatsrcdir, path)
        if (path.startswith('.')): continue
        if (path.startswith('master-darkflat.fits')): continue
        # check if current path is a file
        if os.path.isfile(dfile):
            darkflatfiles.append(path)
    darkflatfiles.sort()
# Do flats, if any
flatfiles = []
if flatsrcdir:
    for path in os.listdir(flatsrcdir):
        dfile = os.path.join(flatsrcdir, path)
        if (path.startswith('.')): continue
        if (path.startswith('master-flat.fits')): continue
        # check if current path is a file
        if os.path.isfile(dfile):
            flatfiles.append(path)
    flatfiles.sort()

# Build dark flat frame, if we have any to work with
darkflat = fits.HDUList()
if len(darkflatfiles) > 0:
    logger.info(f"Processing {len(darkflatfiles)} dark-flats")
    darkflat = buildMedianStack(darkflatsrcdir, darkflatfiles, "master-darkflat.fits")
    calstat += "B"

# Build dark frame, if we have any to work with
dark = fits.HDUList()
if len(darkfiles) > 0:
    logger.info(f"Processing {len(darkfiles)} darks")
    dark = buildMedianStack(darksrcdir, darkfiles, "master-dark.fits")
    calstat += "D"

# Build flat frame, if we have any to work with
flat = fits.HDUList()
if len(flatfiles) > 0:
    logger.info(f"Processing {len(flatfiles)} flats")
    normflataccum = buildMasterFlatStack(flatsrcdir, flatfiles, "master-flat.fits", darkflat)
    calstat += "F"
    
printfirst = True
cnt = 0

for f in lightfiles:
    try:
        lfile = os.path.join(sciencesrcdir, f)
        newfname = "science-{0:05d}.fits".format(cnt)
        newfits = os.path.join(sciencepath, newfname)
        # Load file into list of HDU list 
        with fits.open(lfile) as hduList:
            img_dtype = hduList[0].data.dtype
            # First, calibrate image
            if len(dark) > 0:
                # Clamp the data with the dark from below, so we can subtract without rollover
                np.maximum(hduList[0].data, dark[0].data, out=hduList[0].data)
                # And subtract the dark
                np.subtract(hduList[0].data, dark[0].data, out=hduList[0].data)
            # If we have flat, apply it
            if len(flat) > 0:
                hduList[0].data = hduList[0].data.astype(np.float32) / normflataccum
            hduList[0].data = hduList[0].data.astype(img_dtype)
            if calstat != "":
                hdrList[0].header.set("CALSTAT", calstat)
            hduList.writeto(newfits, overwrite=True)
            cnt = cnt + 1
    except OSError as e:
        logger.error("Error: file %s - %s (%s)" % (f, e.__class__, e))     

logger.info("Processed %d out of %d files into destination '%s'" % (cnt, len(lightfiles), outputdir))

