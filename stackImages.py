import argparse
from genericpath import isfile
import os
import pathlib
import shutil
from tkinter import E
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
import cv2
import numpy as np
import subprocess
import warnings
# UTC to BJD converter import
from barycorrpy import utc_tdb
from pandas import isna

warnings.simplefilter('ignore', category=FITSFixedWarning)

def runsolving(ra, dec, infile, outfile):
    try:
        rslt = subprocess.run(["solve-field", infile,
            "--no-plots", "--overwrite",
            "--ra", str(ra),
            "--dec", str(dec),
            "--radius", "5",
            "--fits-image", "--guess-scale",
            "--new-fits", outfile ], 
            timeout=30, capture_output=True)
        if rslt.returncode != 0:
            print("Error solving %s - skipping" % f)
            return False
        return True
    except subprocess.TimeoutExpired:
        print("Timeout solving %s - skipping" % f)
        return False

def runstacking(ra, dec, fitsfiles, stackfile, mjdobs, mjdend): 
    print("Stacking %d images into %s" % (len(fitsfiles), stackfile))
    tmpst = os.path.join(tmppath, "tmpstack.fits")
    stackargs = [ "SWarp", 
        "-IMAGEOUT_NAME", tmpst, 
        "-WRITE_XML", "N",
        "-RESAMPLE_DIR", tmppath,
        "-COPY_KEYWORDS", "OBJECT,ORIGIN,MINSYET,TELESCOP,INSTUME,SERIALNB,TIMEUNIT,LATITUDE,LONGITUD,GAIN,GAINDB,ALTITUDE,CMOSTEMP,OBSMODE,DATE,SOFTVER" ]                              
    stackargs.extend(fitsfiles)
    rslt = subprocess.run(stackargs, capture_output=True)
    if rslt.returncode != 0:
        print("Error stacking %s" % stackfile)
        return False
    # And resolve new file
    if runsolving(ra, dec, tmpst, stackfile):
        # And add MJD fields for stack
        with fits.open(stackfile) as hduList:
            hduList[0].header['MJD-OBS'] = mjdobs
            hduList[0].header['MJD-END'] = mjdend
            hduList[0].header['MJD-MID'] = (mjdobs + mjdend) / 2
            hduList.writeto(stackfile, overwrite=True)
            print("Stacked file %s" % stackfile)
        return True
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
parser.add_argument("-i", "--input", help = "Source image directory", required=True);
# Adding output argument
parser.add_argument("-o", "--output", help = "Output stacked image directory") 
parser.add_argument("--stacktime", help = "Number of mintutes to stack (default 2)") 

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
stacktime = 2
if args.stacktime:
    stacktime = int(args.stacktime)
    
# Make output directory, if needed
pathlib.Path(outputdir).mkdir(parents=True, exist_ok=True)
tmppath = os.path.join(outputdir, "tmp")
pathlib.Path(tmppath).mkdir(parents=True, exist_ok=True)

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
            tobs = hduList[0].header['MJD-OBS'] * 24 * 60 # MJD in minutes
            # End of accumulator?
            if (len(timeaccumlist) > 0) and ((timeaccumstart + stacktime) < tobs):
                stackout = os.path.join(outputdir, "stack-%04d.fits" % stackedcnt)
                tmpout = os.path.join(tmppath, "tmpstack.fits")
                try: 
                    runstacking(timeaccumra, timeaccumdec, timeaccumlist, stackout, mjdobs, mjdend)
                    stackedcnt = stackedcnt + 1
                except OSError as e:
                    print("Error: stacking file %s - %s (%s)" % (stackout, e.__class__, e))   
                timeaccumlist = []
                timeaccumstart = 0
            else:
                # If first one to accumulate, save start time and RA/Dec
                if (len(timeaccumlist) == 0):
                    timeaccumstart = tobs
                    timeaccumra = hduList[0].header['FOVRA']
                    timeaccumdec = hduList[0].header['FOVDEC']
                    mjdobs = hduList[0].header['MJD-OBS']
                timeaccumlist.append(lfile)
                mjdend = hduList[0].header['MJD-END']
        cnt = cnt + 1
    except OSError as e:
        print("Error: file %s - %s (%s)" % (f, e.__class__, e))     
# Final accumulator?
if len(timeaccumlist) > 0:
    stackout = os.path.join(outputdir, "stack-%04d.fits" % stackedcnt)
    try: 
        runstacking(timeaccumra, timeaccumdec, timeaccumlist, stackout, mjdobs, mjdend)
        stackedcnt = stackedcnt + 1
    except OSError as e:
        print("Error: stacking file %s - %s (%s)" % (stackout, e.__class__, e))     

print("Processed %d images into %d stacked images" % (cnt, stackedcnt))

