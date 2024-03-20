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
parser.add_argument("-i", "--input", help = "Input directory");
# Read arguments from command line
try:
    args = parser.parse_args()
except argparse.ArgumentError:
    os.exit(1)

inputsrcdir='input'
if args.input:
    inputsrcdir = args.input 
# Go through the lights
lightfiles = []
for path in os.listdir(inputsrcdir):
    if (path.startswith('.')): continue
    if not path.lower().endswith('.fit'):
        continue
    dfile = os.path.join(inputsrcdir, path)
    # check if current path is a file
    if os.path.isfile(dfile):
        lightfiles.append(path)
lightfiles.sort()

for f in lightfiles:
    try:
        lfile = os.path.join(inputsrcdir, f)
        # Load file into list of HDU list 
        with fits.open(lfile) as hduList:
            framebuf = hduList[0].data
            xdim, ydim = framebuf.shape
            xmin = int(0.4 * xdim)
            ymin = int(0.4 * ydim)
            xmax = int(0.6 * xdim)
            ymax = int(0.6 * ydim)
            middlebuf = framebuf[xmin:xmax, ymin:ymax]
            print(f"{f}: {xdim},{ydim} mean={np.mean(middlebuf)}")

    except OSError as e:
        logger.error("Error: file %s - %s (%s)" % (f, e.__class__, e))     


