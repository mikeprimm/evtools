import argparse
from genericpath import isfile
import os
import pathlib
import logging
import time
from astropy.io import fits
from astropy.wcs import FITSFixedWarning
from astropy.time import Time
import astroalign as aa
import numpy as np
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


def scaleUp(ary: np.array, scale: int):
   if scale == 1:
       return np.copy(ary)
   return ary.repeat(scale, axis=0).repeat(scale, axis=1)

def scaleDown(ary: np.array, scale: int, dtype: np.dtype):
    if scale == 1:
       return ary.astype(dtype)
    rslt = np.zeros([ ary.shape[0] // scale, ary.shape[1] // scale ], dtype=dtype)
    for x in range(scale):
        for y in range(scale):
            rslt += ary[x::scale,y::scale]
    rslt /= (scale * scale)
    return rslt

# create logger
logger = logging.getLogger('processExoplanetData')
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
# Add flags (default is green)
parser.add_argument('-r', "--red", action='store_true')
parser.add_argument('-g', "--green", action='store_true')
parser.add_argument('-b', "--blue", action='store_true')
parser.add_argument("-G", "--gray", action='store_true')
parser.add_argument('-B', "--bin", action='store_true')
parser.add_argument("-st", "--stacktime", help = "Number of seconds to stack (default 120)") 
parser.add_argument("-sm", "--stackmin", help = "Minimum number of frames per stack (default 5)") 
parser.add_argument('-ss', "--supersample", help = "Supersample scale")
parser.add_argument('-sk', "--skip", help = "Skip frames (every Nth frame)")

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
ch2 = logging.FileHandler(os.path.join(outputdir, 'processExoplanetData.log'), encoding='utf-8', mode='w')
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
if args.darkflats:
   darkflatpath = os.path.join(outputdir, "darkflats")
   pathlib.Path(darkflatpath).mkdir(parents=True, exist_ok=True)
if args.flats:
   flatpath = os.path.join(outputdir, "flats")
   pathlib.Path(flatpath).mkdir(parents=True, exist_ok=True)
stacktime = 120
if args.stacktime:
    stacktime = int(args.stacktime)
supersample = 2
if args.supersample:
    supersample = int(args.supersample) 
stackmin = 5
if args.stackmin:
    stackmin = int(args.stackmin) 
skipcnt = 1
if args.skip:
    skipcnt = int(args.skip) 

doRed = False
doGreen = False
doBlue = False
doGray = False
doBin = False
if args.red:
    doRed = True
    print("Produce red channel FITS files")
elif args.gray:
    doGray = True
    print("Produce grayscale channel FITS files")
elif args.blue:
    doBlue = True
    print("Produce blue channel FITS files")
else:
    doGreen = True
    print("Produce green channel FITS files")
if args.bin:
    doBin = True
    print("2x2 bin files (1/2 resolution)")

darkfiles = []
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

# Build dark frame, if we have any to work with
dark = fits.HDUList()
if len(darkfiles) > 0:
    logger.info(f"Processing {len(darkfiles)} darks")
    dark = buildMedianStack(darksrcdir, darkfiles, "master-dark.fits")
# Build dark flat frame, if we have any to work with
darkflat = fits.HDUList()
if len(darkflatfiles) > 0:
    logger.info(f"Processing {len(darkflatfiles)} dark-flats")
    darkflat = buildMedianStack(darkflatsrcdir, darkflatfiles, "master-darkflat.fits")

# Build flat frame, if we have any to work with
flat = fits.HDUList()
if len(flatfiles) > 0:
    logger.info(f"Processing {len(flatfiles)} flats")
    normflataccum = buildMasterFlatStack(flatsrcdir, flatfiles, "master-flat.fits", darkflat)
    
cnt = 0
stackedcnt = 0
timeaccumcnt = 0
timeaccumstart = 0
timeaccumlast = 0
firstframe = None
accumulatorframe = None
accumulatorcounts = None
accumfname = None
mjdobs = 0
mjdend = 0

lastidx = len(lightfiles) - 1

bayerpat = None

for idx in range(lastidx + 1):
    if (skipcnt > 0) and ((idx % skipcnt) != 0):
        continue
    f = lightfiles[idx]
    try:
        lfile = os.path.join(sciencesrcdir, f)
        # Load file into list of HDU list 
        with fits.open(lfile) as hduList:
            data = hduList[0].data.astype(np.float64)
            if bayerpat is None:
                if 'BAYERPAT' in hduList[0].header:
                    bayerpat = hduList[0].header['BAYERPAT']
                    bayerpat = bayerpat[0:4]
                else:
                    bayerpat = "RGGB"
                print(f"Bayer pattern={bayerpat}")
            # First, calibrate image
            if len(dark) > 0:
                # Clamp the data with the dark from below, so we can subtract without rollover
                np.maximum(data, dark[0].data, out=data)
                # And subtract the dark
                np.subtract(data, dark[0].data, out=data)
            # If we have flat, apply it
            if len(flat) > 0:
                data /= normflataccum
            # Process color channel
            # If making red, split out red channel
            if doBin:
                hduList[0].header.set('FOVXREF', hduList[0].header['FOVXREF'] // 2)
                hduList[0].header.set('FOVYREF', hduList[0].header['FOVYREF'] // 2)
                if doRed:
                    data = data[::2,::2]
                # If making green, split out green channels
                elif doGreen:
                    data = (data[1::2,::2] + data[::2,1::2]) / 2
                # If making blue, split out blue channels
                elif doBlue:
                    data = data[1::2,1::2]
                else:
                    print("gray not supported with binning")                    
            else:
                data = demosaicing_CFA_Bayer_bilinear(data, bayerpat)
                if doRed:
                    data = data @ np.array([ 1, 0, 0 ])
                # If making green, split out green channels
                elif doGreen:
                    data = data @ np.array([ 0, 1, 0 ])
                # If making blue, split out blue channels
                elif doBlue:
                    data = data @ np.array([ 0, 0, 1 ])
                # If making grayscale
                elif doGray:
                    data = data @ np.array([ 0.2125, 0.7154, 0.0721 ]);
            hduList[0].header.remove('BAYERPAT')
            # And stack image
            mjdstart = hduList[0].header['MJD-OBS']
            mjdend = hduList[0].header['MJD-END']
            dateend = hduList[0].header['DATE-END']
            # If past end of our limit, or last image, flush
            tobs = mjdend * 24 * 60 * 60  # MJD in seconds
            if (timeaccumcnt > 0) and (((timeaccumstart + stacktime) < tobs) or (idx == lastidx)):
                if timeaccumcnt < stackmin:
                    print(f"Skip frames from {mjdobs} to {mjdend}")
                else:
                    print(f"Accumulated {timeaccumcnt} frames from {mjdobs} to {mjdend} into frame {stackedcnt}")
                    with fits.open(accumfname) as hduStackList:
                        hduStackList[0].header.set("MJD-OBS", mjdobs)
                        hduStackList[0].header.set("MJD-MID", (mjdobs + mjdend) / 2)
                        hduStackList[0].header.set("MJD-END", mjdend)
                        hduStackList[0].header.set("EXPTIME", (mjdend - mjdobs) * 24 * 3600)
                        hduStackList[0].header.set("DATE-OBS", datestart)
                        hduStackList[0].header.set("DATE-END", dateend)
                        hduStackList[0].header.set("BZERO", 0)
                        hduStackList[0].header.set("BSCALE", 1)
                        stime = Time(datestart)
                        etime = Time(dateend)
                        mtime = Time((stime.jd + etime.jd) / 2, format="jd", scale="tt")
                        mtime.format = "isot"
                        hduStackList[0].header.set("DATE-AVG", mtime.to_string())
                        accumulatorframe /= accumulatorcounts
                        hduStackList[0].data = scaleDown(accumulatorframe, supersample, np.float32)
                        hduStackList.writeto(os.path.join(outputdir, "stack-{0:05d}.fits".format(stackedcnt)), overwrite=True)
                    stackedcnt = stackedcnt + 1
                # Reset accumulator
                timeaccumcnt = 0
                timeaccumstart = 0
            # If first file of new accumulator, add it
            if timeaccumcnt == 0:
                accumfname = lfile
                mjdobs = mjdstart
                datestart = hduList[0].header['DATE-OBS']
                timeaccumstart = mjdstart * 24 * 60 * 60
                firstframe = scaleUp(data, supersample)
                accumulatorframe = np.copy(firstframe)
                accumulatorcounts = np.ones(accumulatorframe.shape, dtype=np.int32)
                timeaccumcnt += 1
            else:
                # Find transformed image to align with first one
                try:
                    registered_image, footprint = aa.register(data, firstframe)
                    accumulatorcounts[footprint == False] += 1
                    accumulatorframe[footprint == False] += registered_image[footprint == False]
                    timeaccumcnt += 1
                except (ValueError, aa.MaxIterError, IndexError) as e:
                    print("Error: Cannot find transform for file %s - %s (%s)" % (f, e.__class__, e))                         
        cnt = cnt + 1
    except OSError as e:
        print("Error: file %s - %s (%s)" % (f, e.__class__, e))     

logger.info("Processed %d out of %d files into destination '%s'" % (cnt, len(lightfiles), outputdir))

