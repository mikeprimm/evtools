import argparse
from genericpath import isfile
import os
import pathlib
from tkinter import E
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.coordinates import SkyCoord
import cv2
import numpy as np
import subprocess
import warnings

from pandas import isna

warnings.simplefilter('ignore', category=FITSFixedWarning)

def runsolving(ra, dec, infile, outfile):
    try:
        rslt = subprocess.run(["solve-field", infile,
            "--no-plots", "--overwrite",
            "--ra", str(ra),
            "--dec", str(dec),
            "--radius", "5",
            "--new-fits", outfile ], 
            timeout=30, capture_output=True)
        if rslt.returncode != 0:
            print("Error solving %s - skipping" % f)
            return False
        return True
    except subprocess.TimeoutExpired:
        print("Timeout solving %s - skipping" % f)
        return False

# Initialize parser
parser = argparse.ArgumentParser()
# Add input argument
parser.add_argument("-d", "--darks", help = "Dark files source directory");
parser.add_argument("-s", "--science", help = "Science files source directory");
# Adding output argument
parser.add_argument("-o", "--output", help = "Output directory") 
# Add flags (default is grey - others for monochrome)
parser.add_argument('-r', "--red", action='store_true')
parser.add_argument('-g', "--green", action='store_true')
parser.add_argument('-b', "--blue", action='store_true')
parser.add_argument("-G", "--gray", action='store_true')

# Read arguments from command line
try:
    args = parser.parse_args()
except argparse.ArgumentError:
    os.exit(1)
outputdir='output'
if args.output: 
    outputdir = args.output
darksrcdir='darks'
if args.darks:
    darksrcdir = args.darks 
sciencesrcdir='science'
if args.science:
    sciencesrcdir = args.science
# Make output directory, if needed
pathlib.Path(outputdir).mkdir(parents=True, exist_ok=True)
darkpath = os.path.join(outputdir, "darks")
pathlib.Path(darkpath).mkdir(parents=True, exist_ok=True)
sciencepath = os.path.join(outputdir, "science")
pathlib.Path(sciencepath).mkdir(parents=True, exist_ok=True)
tmppath = os.path.join(outputdir, "tmp")
pathlib.Path(tmppath).mkdir(parents=True, exist_ok=True)

togray = False
coloridx = 0
if args.red:
    coloridx = 2   # Red
    print("Produce red channel FITS files")
elif args.green:
    coloridx = 1   # Green
    print("Produce green channel FITS files")
elif args.blue:
    coloridx = 0   # Blue
    print("Produce blue channel FITS files")
else:
    togray = True
    print("Produce grayscale FITS files")
darkfiles = []
# Go through the darks
for path in os.listdir(darksrcdir):
    dfile = os.path.join(darksrcdir, path)
    # check if current path is a file
    if os.path.isfile(dfile):
        darkfiles.append(path)
darkfiles.sort()
# Go through the lights
lightfiles = []
for path in os.listdir(sciencesrcdir):
    dfile = os.path.join(sciencesrcdir, path)
    # check if current path is a file
    if os.path.isfile(dfile):
        lightfiles.append(path)
lightfiles.sort()
dark = fits.HDUList()

# Build dark frame, if we have any to work with
if len(darkfiles) > 0:
    for f in darkfiles:
        try:
            dfile = os.path.join(darksrcdir, f)
            # Load file into list of HDU list 
            with fits.open(dfile) as hduList:
                # Use first one as base
                if len(dark) == 0:
                    darkaccum = np.zeros((0,) + hduList[0].data.shape)
                    dark.append(hduList[0].copy())
                darkaccum = np.append(darkaccum, [ hduList[0].data ], axis=0)
                hduList.writeto(os.path.join(darkpath, f), overwrite=True)
        except OSError:
            print("Error: file %s" % f)        
    # Now compute median for each pixel
    darkaccum = np.median(darkaccum, axis=0)
    dark[0].data = darkaccum.astype(np.uint16)
    # And write output dark
    dark.writeto(os.path.join(darkpath, "master-dark.fits"), overwrite=True)

cnt = 0
solvedcnt = 0
timeaccumlist = []
timeaccumstart = 0
timeaccumra = 0
timeaccumdec = 0
mjdobs = 0
mjdend = 0
stackedcnt = 0

for f in lightfiles:
    try:
        lfile = os.path.join(sciencesrcdir, f)
        # Load file into list of HDU list 
        with fits.open(lfile) as hduList:
            # First science? get center as target
            if (cnt == 0):
                fov = SkyCoord(hduList[0].header['FOVRA'], hduList[0].header['FOVDEC'], frame='icrs', unit='deg')
                print("center of FOV for first science: RA={0}, DEC={1}".format(hduList[0].header['FOVRA'], hduList[0].header['FOVDEC']))
            # First, calibrate image
            if len(dark) > 0:
                # Clamp the data with the dark from below, so we can subtract without rollover
                np.maximum(hduList[0].data, dark[0].data, out=hduList[0].data)
                # And subtract the dark
                np.subtract(hduList[0].data, dark[0].data, out=hduList[0].data)
            # Now debayer into grayscale                
            if togray:
                dst = cv2.cvtColor(hduList[0].data, cv2.COLOR_BayerRG2GRAY)
                for idx, val in enumerate(dst):
                    hduList[0].data[idx] = val
            else:
                # Demosaic the image
                dst = cv2.cvtColor(hduList[0].data, cv2.COLOR_BayerRG2BGR)
                for idx, val in enumerate(dst):
                    hduList[0].data[idx] = val[:,coloridx]
            rslt = True
            newfits = os.path.join(sciencepath, "science-{0:05d}.fits".format(cnt))
            # Write to temporary file so that we can run solve-field to
            # set WCS data
            hduList.writeto(os.path.join(tmppath, "tmp.fits"), overwrite=True)
            # Now run solve-field to generate final file
            rslt = runsolving(hduList[0].header['FOVRA'], hduList[0].header['FOVDEC'],
                os.path.join(tmppath, "tmp.fits"), newfits )
            if rslt == True:
                # Read new file - see if we are still in frame
                with fits.open(newfits) as hduListNew:
                    w = WCS(hduListNew[0].header)
                    shape = hduList[0].data.shape
                    x, y = w.world_to_pixel(fov)
                    print("Solved %s:  target at %f, %f" % (newfits, x, y))
                    # If out of range, drop the frame
                    if (x < 0) or (x >= shape[0]) or (y < 0) or (y >= shape[1]):
                        rslt = False;
                        print("Discarding %s - target out of frame" % fname)                    
            if rslt == False:
                print("Error solving %s - skipping" % f)
            else:
                tobs = hduList[0].header['MJD-OBS'] * 24 * 60 # MJD in minutes
                solvedcnt = solvedcnt + 1
            cnt = cnt + 1
    except OSError as e:
        print("Error: file %s - %s (%s)" % (f, e.__class__, e))     

print("Processed %d out of %d files and %d stacks into destination '%s'" % (solvedcnt, cnt, stackedcnt, outputdir))

