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

def calcAltAz(ra, dec, lat, lon, alt, mjdtime):
    pointing = SkyCoord(str(ra) + " " + str(dec), unit=(u.deg, u.deg), frame='icrs')
    location = EarthLocation.from_geodetic(lat=lat * u.deg, lon=lon * u.deg, height=alt)
    atime = Time(mjdtime, format='mjd', scale='utc', location=location)
    pointingAltAz = pointing.transform_to(AltAz(obstime=atime, location=location))
    return pointingAltAz

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("--ra", help = "Target RA", required = True)
parser.add_argument("--dec", help = "Target Dec", required = True)
parser.add_argument("--obslat", help = "Observatory Latitude (deg)", required = True)
parser.add_argument("--obslon", help = "Observatory Longitude (deg east)", required = True)
parser.add_argument("--obsalt", help = "Observatory Altitude (meters)", required = True)
parser.add_argument("--bjdtdb", help = "Time (BJD_TDB)", required = True)
# Read arguments from command line
try:
    args = parser.parse_args()
except argparse.ArgumentError:
    os.exit(1)

logger.info("args=%s" % args)

obsAltitude = None
obsLatitude = None
obsLongitude = None
bjdtdb = None
if args.obslat:
    obsLatitude = float(args.obslat)
if args.obslon:
    obsLongitude = float(args.obslon)
if args.obsalt:
    obsAltitude = float(args.obsalt)
if args.bjdtdb:
    bjdtdb = float(args.bjdtdb)
# Get target
if args.dec:
    target = SkyCoord(args.ra, args.dec, frame='icrs', unit=(u.hourangle, u.deg))
    targetRA = target.ra.deg
    targetDec = target.dec.deg
    logger.info("Target coords: RA=%d:%d:%f, Dec=%s%d:%d:%f" % (target.ra.hms.h, target.ra.hms.m, target.ra.hms.s, '+' if       target.dec.signed_dms.sign >= 0 else '-', target.dec.signed_dms.d, target.dec.signed_dms.m, target.dec.signed_dms.s))
else:
    targetRA = None
    targetDec = None

# Add alt-az for target
altaz = calcAltAz(targetRA, targetDec, obsLatitude, obsLongitude, obsAltitude, bjdtdb)
logger.info(f"alt-az={altaz}");
logger.info(f"airmass={float(altaz.secz)}");

