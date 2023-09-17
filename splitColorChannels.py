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
logger = logging.getLogger('splitColorChannels')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

# Initialize parser
parser = argparse.ArgumentParser()
# Add input argument
parser.add_argument("-s", "--science", help = "Science files source directory");
# Adding output argument
parser.add_argument("-o", "--output", help = "Output directory") 
# Add flags (default is grey - others for monochrome)
parser.add_argument('-r', "--red", action='store_true')
parser.add_argument('-g', "--green", action='store_true')
parser.add_argument('-b', "--blue", action='store_true')
parser.add_argument('-A', "--all", action='store_true')
parser.add_argument('-B', "--bin", action='store_true')
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
ch2 = logging.FileHandler(os.path.join(outputdir, 'splitColorChannels.log'), encoding='utf-8', mode='w')
ch2.setLevel(logging.DEBUG)
ch2.setFormatter(formatter)
logger.addHandler(ch2)

logger.info("args=%s" % args)

sciencesrcdir='science'
if args.science:
    sciencesrcdir = args.science

# Make output directory, if needed
pathlib.Path(outputdir).mkdir(parents=True, exist_ok=True)

doRed = False
doGreen = False
doBlue = False
doBin = False
if args.red:
    doRed = True
    print("Produce red channel FITS files")
elif args.all:
    doRed = True
    doGreen = True
    doBlue = True
    print("Produce red, green, and blue channel FITS files")
elif args.blue:
    print("Produce blue channel FITS files")
else:
    doGreen = True
    print("Produce green channel FITS files")
if args.bin:
    doBin = True
    print("2x2 bin files (1/2 resolution)")

# Go through the lights
lightfiles = []
for path in os.listdir(sciencesrcdir):
    if (path.startswith('.')): continue
    dfile = os.path.join(sciencesrcdir, path)
    # check if current path is a file
    if os.path.isfile(dfile):
        lightfiles.append(path)
lightfiles.sort()
    
cnt = 0

for f in lightfiles:
    try:
        lfile = os.path.join(sciencesrcdir, f)
        # Load file into list of HDU list 
        with fits.open(lfile) as hduList:
            data = hduList[0].data
            img_dtype = hduList[0].data.dtype
            # If making red, split out red channel
            if doBin:
                hduList[0].header.set('FOVXREF', hduList[0].header['FOVXREF'] // 2)
                hduList[0].header.set('FOVYREF', hduList[0].header['FOVYREF'] // 2)
                if doRed:
                    hduList[0].data = data[::2,::2]
                    hduList.writeto(os.path.join(outputdir, "red-{0:05d}.fits".format(cnt)), overwrite=True)
                # If making green, split out green channels
                if doGreen:
                    hduList[0].data = (data[1::2,::2] + data[::2,1::2]) / 2
                    hduList[0].data = hduList[0].data.astype(img_dtype)
                    hduList.writeto(os.path.join(outputdir, "green-{0:05d}.fits".format(cnt)), overwrite=True)
                # If making blue, split out blue channels
                if doBlue:
                    hduList[0].data = data[1::2,1::2]
                    hduList.writeto(os.path.join(outputdir, "blue-{0:05d}.fits".format(cnt)), overwrite=True)
            else:
                new_image_data = demosaicing_CFA_Bayer_bilinear(hduList[0].data, "RGGB")
                if doRed:
                    hduList[0].data = new_image_data @ np.array([ 1, 0, 0 ])
                    hduList[0].data = hduList[0].data.astype(img_dtype)
                    hduList.writeto(os.path.join(outputdir, "red-{0:05d}.fits".format(cnt)), overwrite=True)
                # If making green, split out green channels
                if doGreen:
                    hduList[0].data = new_image_data @ np.array([ 0, 1, 0 ])
                    hduList[0].data = hduList[0].data.astype(img_dtype)
                    hduList.writeto(os.path.join(outputdir, "green-{0:05d}.fits".format(cnt)), overwrite=True)
                # If making blue, split out blue channels
                if doBlue:
                    hduList[0].data = new_image_data @ np.array([ 0, 0, 1 ])
                    hduList[0].data = hduList[0].data.astype(img_dtype)
                    hduList.writeto(os.path.join(outputdir, "blue-{0:05d}.fits".format(cnt)), overwrite=True)

            cnt = cnt + 1
    except OSError as e:
        logger.error("Error: file %s - %s (%s)" % (f, e.__class__, e))     

logger.info("Processed %d out of %d files into destination '%s'" % (cnt, len(lightfiles), outputdir))

