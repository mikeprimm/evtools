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
parser.add_argument("-d", "--darks", help = "Dark files source directory");
parser.add_argument("-s", "--science", help = "Science files source directory");
parser.add_argument("-df", "--darkflats", help = "Dark flat files source directory");
parser.add_argument("-f", "--flats", help = "Flat files source directory");
# Adding output argument
parser.add_argument("-o", "--output", help = "Output directory") 
# Add flags (default is grey - others for monochrome)
parser.add_argument('-r', "--red", action='store_true')
parser.add_argument('-g', "--green", action='store_true')
parser.add_argument('-b', "--blue", action='store_true')
parser.add_argument("-G", "--gray", action='store_true')
parser.add_argument("-B", "--blueblock", action='store_true')
parser.add_argument('-A', "--all", action='store_true')
parser.add_argument("-N", "--nosolve", action='store_true')
parser.add_argument("--target", help = "Target Name (EXOFOP)")
parser.add_argument("--ra", help = "Target RA")
parser.add_argument("--dec", help = "Target Dec")
parser.add_argument("--c1ra", help = "Comparison 1 RA")
parser.add_argument("--c1dec", help = "Comparison 1 Dec")
parser.add_argument("--obslat", help = "Observatory Latitude (deg)")
parser.add_argument("--obslon", help = "Observatory Longitude (deg east)")
parser.add_argument("--obsalt", help = "Observatory Altitude (meters)")
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
ch2 = logging.FileHandler(os.path.join(outputdir, 'processExoplanetTransit.log'), encoding='utf-8', mode='w')
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

borderbuffer = 20
if args.borderbuffer:
    borderbuffer = args.borderbuffer
obsAltitude = None
obsLatitude = None
obsLongitude = None
if args.obslat:
    obsLatitude = float(args.obslat)
if args.obslon:
    obsLongitude = float(args.obslon)
if args.obsalt:
    obsAltitude = float(args.obsalt)

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

c1 = None
if args.c1ra:
   c1 = SkyCoord(args.c1ra, args.c1dec, frame='icrs', unit=(u.hourangle, u.deg))
   c1RA = c1.ra.deg
   c1Dec = c1.dec.deg
   logger.info("Comparison 1 coords: RA=%d:%d:%f, Dec=%s%d:%d:%f" % (c1.ra.hms.h, c1.ra.hms.m, c1.ra.hms.s, '+' if c1.dec.signed_dms.sign >= 0 else '-', c1.dec.signed_dms.d, c1.dec.signed_dms.m, c1.dec.signed_dms.s))

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

togray = False
blueblock = False
tobayer = False
toall = False
coloridx = 0
fltname='bayer'
solve = True

if args.red:
    # Red
    colorweight = np.array([ 1, 0, 0 ])
    fltname='CR'
    print("Produce red channel FITS files")
elif args.green:
    # Green
    colorweight = np.array([ 0, 1, 0 ])
    fltname='TG'
    print("Produce green channel FITS files")
elif args.blue:
    # Blue
    colorweight = np.array([ 0, 0, 1 ]);
    fltname='TB'
    print("Produce blue channel FITS files")
elif args.blueblock:
    colorweight = np.array([ 0.2125, 0.7154, 0 ]);
    fltname='CBB'
    print("Produce blueblocked grayscale FITS files")
elif args.gray:
    colorweight = np.array([ 0.2125, 0.7154, 0.0721 ]);
    fltname='CV'
    print("Produce grayscale FITS files")
else:
    tobayer = True
    fltname='BAYER'
    if args.all:
        toall = True
        print("Produce red, green, blue, gray channel FITS files")
        redpath = os.path.join(sciencepath, "red")
        greenpath = os.path.join(sciencepath, "green")
        bluepath = os.path.join(sciencepath, "blue")
        graypath = os.path.join(sciencepath, "gray")
        pathlib.Path(redpath).mkdir(parents=True, exist_ok=True)
        pathlib.Path(greenpath).mkdir(parents=True, exist_ok=True)
        pathlib.Path(bluepath).mkdir(parents=True, exist_ok=True)
        pathlib.Path(graypath).mkdir(parents=True, exist_ok=True)
    else:
        logger.info("Produce Bayer FITS files")

if args.nosolve:
    solve = False
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

cnt = 0
timeaccumlist = []
timeaccumstart = 0
timeaccumra = 0
timeaccumdec = 0
mjdobs = 0
mjdend = 0
printfirst = True

for f in lightfiles:
    try:
        lfile = os.path.join(sciencesrcdir, f)
        newfname = "science-{1}-{0:05d}.fits".format(cnt, fltname)
        newfits = os.path.join(sciencepath, newfname)
        # Load file into list of HDU list 
        with fits.open(lfile) as hduList:
            if obsAltitude is None:
                obsAltitude = hduList[0].header['ALTITUDE']
            if obsLatitude is None:
                obsLatitude = hduList[0].header['LATITUDE']
            if obsLongitude is None:
                obsLongitude = hduList[0].header['LONGITUD']
            if cnt == 0:
                logger.info("Observatory: Lat={0} deg, Lon={1} deg, Alt={2} meters".format(obsLatitude, obsLongitude, obsAltitude))
                if target:
                    targetNow = target.apply_space_motion(new_obstime = Time(hduList[0].header['MJD-MID'], format='mjd'))
                    logger.info("Target coords (obs date): RA=%d:%d:%f, Dec=%s%d:%d:%f" % (targetNow.ra.hms.h, targetNow.ra.hms.m, targetNow.ra.hms.s, '+' if targetNow.dec.signed_dms.sign >= 0 else '-', targetNow.dec.signed_dms.d, targetNow.dec.signed_dms.m, targetNow.dec.signed_dms.s))
                    targetRA = targetNow.ra.deg
                    targetDec = targetNow.dec.deg
            img_dtype = hduList[0].data.dtype
            # First, calibrate image
            if len(dark) > 0:
                # Clamp the data with the dark from below, so we can subtract without rollover
                np.maximum(hduList[0].data, dark[0].data, out=hduList[0].data)
                # And subtract the dark
                np.subtract(hduList[0].data, dark[0].data, out=hduList[0].data)
            if solve:
                hduList.writeto(os.path.join(tmppath, "tmp.fits"), overwrite=True)
            else:
                hduList.writeto(newfits, overwrite=True)
        rslt = True
        if solve:
            # Now run solve-field to generate final file
            rslt = runsolving(hduList[0].header['FOVRA'], hduList[0].header['FOVDEC'],
                os.path.join(tmppath, "tmp.fits"), newfits )
        if rslt == False:
            print("Error solving %s - skipping" % f)
            hduList.writeto(os.path.join(badsciencepath, f), overwrite=True)
            continue
        with fits.open(newfits) as hduList:
            # Remember base type
            img_dtype = hduList[0].data.dtype
            # If we have flat, apply it
            if len(flat) > 0:
                hduList[0].data = hduList[0].data.astype(np.float32) / normflataccum
            # Now debayer into grayscale                
            if not tobayer:
                # Demosaic the image
                new_image_data = demosaicing_CFA_Bayer_bilinear(hduList[0].data, "RGGB")
                hduList[0].data = new_image_data @ colorweight
            # Compute BJD times
            if targetRA:
                mjdtimes = np.array([hduList[0].header['MJD-MID']])
                bjdtimes = utc_tdb.JDUTC_to_BJDTDB(mjdtimes + 2400000.5, ra=targetRA, dec=targetDec,
                    lat=obsLatitude, longi=obsLongitude, alt=obsAltitude)[0]
                hduList[0].header.set('BJD_TDB', bjdtimes[0], "barycentric Julian date of the mid obs")
                # Add alt-az for target
                altaz = calcAltAz(targetRA, targetDec, obsLatitude, obsLongitude, obsAltitude, mjdtimes[0])
                # Add AIRMASS
                hduList[0].header.set('AIRMASS', float(altaz.secz))
                # Add RAOBJ2K, DECOBJ2K, SITELAT, SITELONG, ALT_OBJ, AZ_OBJ, ZD_OBJ
                hduList[0].header.set('SITELAT', obsLatitude, "geographic latitude of observatory")
                hduList[0].header.set('SITELONG', obsLongitude, "geographic longitude of observatory")
                hduList[0].header.set('RAOBJ2K', targetRA / 15, "J2000 right ascension of target (hours)")
                hduList[0].header.set('DECOBJ2K', targetDec, "J2000 declination of target (degrees)")
                hduList[0].header.set('ALT_OBJ', float(altaz.alt.deg), "Target altitude at mid-exposure")
                hduList[0].header.set('AZ_OBJ', float(altaz.az.deg), "Target azimuth at mid-exposure")
                hduList[0].header.set('ZD_OBJ', 90.0 - float(altaz.alt.deg), "Target zenith distance at mid-exposure")

            # Add bayer header if leaving as bayer file
            if tobayer:
                hduList[0].header.set('BAYERPAT', 'RGGB')
                hduList[0].header.set('XBAYROFF', 0)
                hduList[0].header.set('YBAYROFF', 0)
            if ('ALTITUDE' in hduList[0].header) == False:
                hduList[0].header.set('ALTITUDE', obsAltitude, "altitude in meters of observing site")
            if ('LATITUDE' in hduList[0].header) == False:
                hduList[0].header.set('LATITUDE', obsLatitude, "latitude in degrees north of observing site")
            if ('LONGITUD' in hduList[0].header) == False:
                hduList[0].header.set('LONGITUD', obsLongitude, "longitude in degrees east of observing site")
            if calstat != "":
                hdrList[0].header.set('CALSTAT', calstat)
            hduList[0].data = hduList[0].data.astype(img_dtype)
            hduList.writeto(newfits, overwrite=True)
        # if solved and we have target ra/dec, check it in field
        if solve and targetRA:
            with fits.open(newfits) as hduList:
                w = WCS(hduList[0].header)
                shape = hduList[0].data.shape
                try:
                    noconv = False
                    x, y = w.world_to_pixel(target)
                    if c1 is None:
                        x1 = x
                        y1 = y
                    else:
                        x1, y1 = w.world_to_pixel(c1)
                except NoConvergence as e:
                    noconv = True
                # If out of range, drop the frame
                if noconv:
                    rslt = False;
                    logger.warning("Rejecting - no convergence on frame %s" % (f))
                    shutil.move(newfits, os.path.join(badsciencepath, f))                    
                elif (x < borderbuffer) or (x >= (shape[0]-borderbuffer)) or (y < borderbuffer) or (y >= (shape[1]-borderbuffer)):
                    rslt = False;
                    logger.warning("Rejecting - target out of frame %s (%f, %f)" % (f, x, y))
                    shutil.move(newfits, os.path.join(badsciencepath, f))
                elif (x1 < borderbuffer) or (x1 >= (shape[0]-borderbuffer)) or (y1 < borderbuffer) or (y1 >= (shape[1]-borderbuffer)):
                    rslt = False;
                    logger.warning("Rejecting - comparison 1 out of frame %s (%f, %f)" % (f, x1, y1))
                    shutil.move(newfits, os.path.join(badsciencepath, f))
                else:
                    if printfirst:
                        logger.info("Pixel coordinate of target: %f, %f" % (x, y))
                        if c1 is not None:
                            logger.info("Pixel coordinate of comparison 1: %f, %f" % (x1, y1))
                        printfirst = False
        # If good result AND we are doing toall, generate different versions
        if rslt and toall:                
            with fits.open(newfits) as hduList:
                hduList[0].header.remove('BAYERPAT')
                hduList[0].header.remove('XBAYROFF')
                hduList[0].header.remove('YBAYROFF')
                new_image_data = demosaicing_CFA_Bayer_bilinear(hduList[0].data, "RGGB")
                # Make red
                red_image_data = new_image_data @ np.array([ 1, 0, 0 ])
                hduList[0].data = red_image_data.astype(img_dtype)
                newfname = "science-{1}-{0:05d}.fits".format(cnt, "CR")
                newfits = os.path.join(redpath, newfname)
                hduList.writeto(newfits, overwrite=True)
                # Make green
                green_image_data = new_image_data @ np.array([ 0, 1, 0 ])
                hduList[0].data = green_image_data.astype(img_dtype)
                newfname = "science-{1}-{0:05d}.fits".format(cnt, "CG")
                newfits = os.path.join(greenpath, newfname)
                hduList.writeto(newfits, overwrite=True)
                # Make blue
                blue_image_data = new_image_data @ np.array([ 0, 0, 1 ])
                hduList[0].data = blue_image_data.astype(img_dtype)
                newfname = "science-{1}-{0:05d}.fits".format(cnt, "CB")
                newfits = os.path.join(bluepath, newfname)
                hduList.writeto(newfits, overwrite=True)
                # Make gray
                gray_image_data = new_image_data @ np.array([ 0.2125, 0.7154, 0.0721 ])
                hduList[0].data = gray_image_data.astype(img_dtype)
                newfname = "science-{1}-{0:05d}.fits".format(cnt, "CV")
                newfits = os.path.join(graypath, newfname)
                hduList.writeto(newfits, overwrite=True)
        if rslt:
            cnt = cnt + 1
    except OSError as e:
        logger.error("Error: file %s - %s (%s)" % (f, e.__class__, e))     

logger.info("Processed %d out of %d files into destination '%s'" % (cnt, len(lightfiles), outputdir))

