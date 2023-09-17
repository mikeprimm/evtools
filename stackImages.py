import argparse
from genericpath import isfile
import os
import pathlib
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.time import Time
import numpy as np
import logging
import astroalign as aa

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
logger = logging.getLogger('stackImages')
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
# Adding output argument
parser.add_argument("-o", "--output", help = "Output stacked image directory") 
parser.add_argument("--stacktime", help = "Number of seconds to stack (default 120)") 

# Read arguments from command line
try:
    args = parser.parse_args()
except argparse.ArgumentError:
    os.exit(1)
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

ch2 = logging.FileHandler(os.path.join(outputdir, 'stackImages.log'), encoding='utf-8', mode='w')
ch2.setLevel(logging.DEBUG)
ch2.setFormatter(formatter)
logger.addHandler(ch2)

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
timeaccumcnt = 0
timeaccumstart = 0
timeaccumlast = 0
firstframe = None
accumulatorframe = None
accumulatorcounts = None

mjdobs = 0
mjdend = 0

lastidx = len(lightfiles) - 1
for idx in range(lastidx + 1):
    f = lightfiles[idx]
    try:
        lfile = os.path.join(inputdir, f)
        # Load file into list of HDU list 
        with fits.open(lfile) as hduList:
            imagedata = hduList[0].data
            mjdstart = hduList[0].header['MJD-OBS']
            if timeaccumcnt == 0:
                mjdobs = mjdstart
                datestart = hduList[0].header['DATE-OBS']
                dateend = hduList[0].header['DATE-END']
                timeaccumstart = mjdstart * 24 * 60 * 60
                firstframe = imagedata
                accumulatorframe = imagedata.astype(np.float64)
                accumulatorcounts = np.ones(accumulatorframe.shape)
                mjdend = hduList[0].header['MJD-END']
                timeaccumcnt += 1
            else:
                # Find transformed image to align with first one
                try:
                    registered_image, footprint = aa.register(imagedata, firstframe)
                    accumulatorcounts[footprint == False] += 1
                    accumulatorframe[footprint == False] += registered_image[footprint == False]
                    mjdend = hduList[0].header['MJD-END']
                    dateend = hduList[0].header['DATE-END']
                    timeaccumcnt += 1
                except ValueError as e:
                    print("Error: Cannot find transform for file %s - %s (%s)" % (f, e.__class__, e))                         
                except aa.MaxIterError as e:
                    print("Error: Cannot find transform for file %s - %s (%s)" % (f, e.__class__, e))                         
            tobs = mjdend * 24 * 60 * 60  # MJD in seconds
            # Past end of accumulator?
            if (timeaccumcnt > 0) and (((timeaccumstart + stacktime) < tobs) or (idx == lastidx)):
                print(f"Accumulated {timeaccumcnt} frames from {mjdobs} to {mjdend} into frame {stackedcnt}")
                hduList[0].header.set("MJD-OBS", mjdobs)
                hduList[0].header.set("MJD-MID", (mjdobs + mjdend) / 2)
                hduList[0].header.set("MJD-END", mjdend)
                hduList[0].header.set("EXPTIME", (mjdend - mjdobs) * 24 * 3600)
                hduList[0].header.set("DATE-OBS", datestart)
                hduList[0].header.set("DATE-END", dateend)
                stime = Time(datestart)
                etime = Time(dateend)
                mtime = Time((stime.jd + etime.jd) / 2, format="jd", scale="tt")
                mtime.format = "isot"
                hduList[0].header.set("DATE-AVG", mtime.to_string())
                hduList[0].data = accumulatorframe / accumulatorcounts
                hduList[0].data = hduList[0].data.astype(np.float32)
                hduList.writeto(os.path.join(outputdir, f"stack-{stackedcnt}.fits"), overwrite=True)
                stackedcnt = stackedcnt + 1
                # Reset accumulator
                timeaccumcnt = 0
                timeaccumstart = 0
        cnt = cnt + 1
    except OSError as e:
        print("Error: file %s - %s (%s)" % (f, e.__class__, e))     

print("Processed %d images into %d stacked images" % (cnt, stackedcnt))

