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

def runsolving(ra, dec, infile, outfile):
    try:
        args = ["solve-field", infile,
            "--no-remove-lines", "--uniformize", "0",
            "--no-plots", "--overwrite",
            "--ra", str(ra),
            "--dec", str(dec),
            "--radius", "5",
            # All Unistellar scopes are between 1 and 2 arcsec per pixel
            "--scale-low", "1",
            "--scale-high", "2",
             "--scale-units", "arcsecperpix",
            "--fits-image",
            "--new-fits", outfile ]
        #logger.info(f"args={args}")
        rslt = subprocess.run(args, 
            timeout=10, capture_output=True)
        if rslt.returncode != 0 or (not isfile(outfile)):
            logger.error("Error solving %s - skipping" % f)
            #logger.error(rslt.stdout)
            #logger.error(rslt.stderr)
            return False
        return True
    except subprocess.TimeoutExpired:
        logger.error("Timeout solving %s - skipping" % f)
        return False

def calcAltAz(ra, dec, lat, lon, alt, mjdtime):
    pointing = SkyCoord(str(ra) + " " + str(dec), unit=(u.deg, u.deg), frame='icrs')
    location = EarthLocation.from_geodetic(lat=lat * u.deg, lon=lon * u.deg, height=alt)
    atime = Time(mjdtime, format='mjd', scale='utc', location=location)
    pointingAltAz = pointing.transform_to(AltAz(obstime=atime, location=location))
    return pointingAltAz

# Initialize parser
parser = argparse.ArgumentParser()
# Add input argument
parser.add_argument("-i", "--input", help = "Input directory");
parser.add_argument("-o", "--output", help = "Output directory") 
parser.add_argument("-t", "--target", help = "Target Name (EXOFOP)")
parser.add_argument("--ra", help = "Target RA")
parser.add_argument("--dec", help = "Target Dec")
parser.add_argument("-bb", "--borderbuffer", help="Border buffer (pixels)", type=int)
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
ch2 = logging.FileHandler(os.path.join(outputdir, 'solveFrames.log'), encoding='utf-8', mode='w')
ch2.setLevel(logging.DEBUG)
ch2.setFormatter(formatter)
logger.addHandler(ch2)

logger.info("args=%s" % args)

inputsrcdir='input'
if args.input:
    inputsrcdir = args.input 
borderbuffer = 20
if args.borderbuffer:
    borderbuffer = args.borderbuffer
# Get target
target = None
targetRA = None
targetDec = None
if args.target:
    tic = exofop_getticid(args.target)
    if tic:
        target, vmag = exofop_getcompositeinfo(tic)
        targetRA = target.ra.deg
        targetDec = target.dec.deg
        logger.info("Target coords (J2000): RA=%d:%d:%f, Dec=%s%d:%d:%f" % (target.ra.hms.h, target.ra.hms.m, target.ra.hms.s, '+' if       target.dec.signed_dms.sign >= 0 else '-', target.dec.signed_dms.d, target.dec.signed_dms.m, target.dec.signed_dms.s))
    else:
        logger.warn("Target not found at EXOFOP: %s" % args.targetcoords)
elif args.dec:
    target = SkyCoord(args.ra, args.dec, frame='icrs', unit=(u.hourangle, u.deg))
    targetRA = target.ra.deg
    targetDec = target.dec.deg
    logger.info("Target coords: RA=%d:%d:%f, Dec=%s%d:%d:%f" % (target.ra.hms.h, target.ra.hms.m, target.ra.hms.s, '+' if       target.dec.signed_dms.sign >= 0 else '-', target.dec.signed_dms.d, target.dec.signed_dms.m, target.dec.signed_dms.s))
# Make output directory, if needed
pathlib.Path(outputdir).mkdir(parents=True, exist_ok=True)
badoutputpath = os.path.join(outputdir, "unsolved")
pathlib.Path(badoutputpath).mkdir(parents=True, exist_ok=True)
tmppath = os.path.join(outputdir, "tmp")
pathlib.Path(tmppath).mkdir(parents=True, exist_ok=True)
# Go through the lights
lightfiles = []
for path in os.listdir(inputsrcdir):
    if (path.startswith('.')): continue
    if not path.lower().endswith('.fits'):
        continue
    dfile = os.path.join(inputsrcdir, path)
    # check if current path is a file
    if os.path.isfile(dfile):
        lightfiles.append(path)
lightfiles.sort()

cnt = 0
for f in lightfiles:
    try:
        lfile = os.path.join(inputsrcdir, f)
        newfname = "science-{0:05d}.fits".format(cnt)
        newfits = os.path.join(outputdir, newfname)
        # Load file into list of HDU list 
        with fits.open(lfile) as hduList:
            if cnt == 0:
                if target:
                    targetNow = target.apply_space_motion(new_obstime = Time(hduList[0].header['MJD-MID'], format='mjd'))
                    logger.info("Target coords (obs date): RA=%d:%d:%f, Dec=%s%d:%d:%f" % (targetNow.ra.hms.h, targetNow.ra.hms.m, targetNow.ra.hms.s, '+' if targetNow.dec.signed_dms.sign >= 0 else '-', targetNow.dec.signed_dms.d, targetNow.dec.signed_dms.m, targetNow.dec.signed_dms.s))
                    targetRA = targetNow.ra.deg
                    targetDec = targetNow.dec.deg
            hduList.writeto(os.path.join(tmppath, "tmp.fits"), overwrite=True)
            # Now run solve-field to generate final file
            rslt = runsolving(hduList[0].header['FOVRA'], hduList[0].header['FOVDEC'],
                os.path.join(tmppath, "tmp.fits"), newfits )
            if rslt == False:
                print("Error solving %s - skipping" % f)
                hduList.writeto(os.path.join(badoutputpath, f), overwrite=True)
                continue                
        with fits.open(newfits) as hduList:
            # Remember base type
            img_dtype = hduList[0].data.dtype
            # Compute BJD times
            if targetRA:
                mjdtimes = np.array([hduList[0].header['MJD-MID']])
                bjdtimes = utc_tdb.JDUTC_to_BJDTDB(mjdtimes + 2400000.5, ra=targetRA, dec=targetDec)[0]
                hduList[0].header.set('BJD_TDB', bjdtimes[0], "barycentric Julian date of the mid obs")
            hduList.writeto(newfits, overwrite=True)
        # if solved and we have target ra/dec, check it in field
        if targetRA:
            with fits.open(newfits) as hduList:
                w = WCS(hduList[0].header)
                shape = hduList[0].data.shape
                try:
                    noconv = False
                    x, y = w.world_to_pixel(target)
                except NoConvergence as e:
                    noconv = True
                # If out of range, drop the frame
                if noconv:
                    rslt = False;
                    logger.warning("Rejecting - no convergence on frame %s" % (f))
                    shutil.move(newfits, os.path.join(badoutputpath, f))                    
                elif (x < borderbuffer) or (x >= (shape[0]-borderbuffer)) or (y < borderbuffer) or (y >= (shape[1]-borderbuffer)):
                    rslt = False;
                    logger.warning("Rejecting - target out of frame %s (%f, %f)" % (f, x, y))
                    shutil.move(newfits, os.path.join(badoutputpath, f))
        if rslt:
            cnt = cnt + 1
    except OSError as e:
        logger.error("Error: file %s - %s (%s)" % (f, e.__class__, e))     

logger.info("Processed %d out of %d files into destination '%s'" % (cnt, len(lightfiles), outputdir))

