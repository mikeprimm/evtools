import argparse
from genericpath import isfile
import os
import pathlib
from astropy.io import fits
import cv2
import csv
import numpy as np
import subprocess

def runsolving(ra, dec, infile, outfile):
    rslt = subprocess.run(["solve-field", infile,
        "--no-plots", "--overwrite",
        "--ra", str(ra),
        "--dec", str(dec),
        "--radius", "5",
        "--new-fits", outfile, "--out", tmppath ], 
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if rslt.returncode != 0:
        print("Error solving %s - skipping" % f)
        return False
    return True

def runstacking(ra, dec, fitsfiles, stackfile): 
    print("Stacking %d images into %s" % (len(fitsfiles), stackfile))
    tmpst = os.path.join(tmppath, "tmpstack.fits")
    stackargs = [ "SWarp", "-IMAGEOUT_NAME", tmpst, "-WRITE_XML", "N" ]
    stackargs.extend(fitsfiles)
    rslt = subprocess.run(stackargs,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if rslt.returncode != 0:
        print("Error stacking %s" % stackfile)
        return False
    # And resolve new file
    return runsolving(ra, dec, tmpst, stackfile)

# Initialize parser
parser = argparse.ArgumentParser()
# Add input argument
parser.add_argument("csvfile", help = "Input CSV file") 
# Adding output argument
parser.add_argument("-o", "--output", help = "Output directory") 
# Add flags (default is grey - others for monochrome)
parser.add_argument('-r', "--red", action='store_true')
parser.add_argument('-g', "--green", action='store_true')
parser.add_argument('-b', "--blue", action='store_true')
parser.add_argument("-G", "--gray", action='store_true')
# Add number of minutes to stack
parser.add_argument("--stacktime", help = "Number of mintutes to stack (default 2)") 

# Read arguments from command line
try:
    args = parser.parse_args()
except argparse.ArgumentError:
    os.exit(1)
outputdir='.'
if args.output: 
    outputdir = args.output
stacktime = 2
if args.stacktime:
    stacktime = int(args.stacktime)
if args.csvfile is None:
    print("CSV filename is required")
    os.exit(1)
if os.path.isfile(args.csvfile) == False:
    print("CSV file not found: ", args.csvfile)
    os.exit(1)
basedir = os.path.dirname(args.csvfile)

# Make output directory, if needed
pathlib.Path(outputdir).mkdir(parents=True, exist_ok=True)
darkpath = os.path.join(outputdir, "darks")
pathlib.Path(darkpath).mkdir(parents=True, exist_ok=True)
sciencepath = os.path.join(outputdir, "science")
pathlib.Path(sciencepath).mkdir(parents=True, exist_ok=True)
tmppath = os.path.join(outputdir, "tmp")
pathlib.Path(tmppath).mkdir(parents=True, exist_ok=True)
stackedpath = os.path.join(outputdir, "stacked")
pathlib.Path(stackedpath).mkdir(parents=True, exist_ok=True)

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
lightfiles = []
with open(args.csvfile, newline='') as csvfile:
    csvreader = csv.DictReader(csvfile)
    for row in csvreader:
        category = row['category']
        if category == 'dark':
            darkfiles.append(row['filename'])
        elif category == 'science':
            lightfiles.append(row['filename']) 
    print("CSV file has %d darks and %d science FITS files" % (len(darkfiles), len(lightfiles)))

dark = fits.HDUList()

# Build dark frame, if we have any to work with
if len(darkfiles) > 0:
    for f in darkfiles:
        try:
            fname = os.path.join(basedir, f)
            # Load file into list of HDU list 
            with fits.open(fname) as hduList:
                # Use first one as base
                if len(dark) == 0:
                    darkaccum = np.zeros(hduList[0].data.shape)
                    dark.append(hduList[0].copy())
                np.add(darkaccum, hduList[0].data, out=darkaccum)
                hduList.writeto(os.path.join(darkpath, f), overwrite=True)
        except OSError:
            print("Error: file %s" % f)        
    # Now compute average for each pixel
    darkaccum = darkaccum // len(darkfiles)
    dark[0].data = darkaccum.astype(np.uint16)
    # And write output dark
    dark.writeto(os.path.join(darkpath, "master-dark.fits"), overwrite=True)

cnt = 0
solvedcnt = 0
timeaccumlist = []
timeaccumstart = 0
timeaccumra = 0
timeaccumdec = 0
stackedcnt = 0
for f in lightfiles:
    try:
        # Load file into list of HDU list 
        fname = os.path.join(basedir, f)
        with fits.open(fname) as hduList:
            print("Processing %s" % fname)
            # First, calibrate image
            if len(dark) > 0:
                # Clamp the data with the dark from below, so we can subtract without rollover
                np.maximum(hduList[0].data, dark[0].data, out=hduList[0].data)
                # And subtract the dark
                np.subtract(hduList[0].data, dark[0].data, out=hduList[0].data)
            # Now debayer into grayscale                
            if togray:
                dst = cv2.cvtColor(hduList[0].data, cv2.COLOR_BayerBG2GRAY)
                for idx, val in enumerate(dst):
                    hduList[0].data[idx] = val
            else:
                # Demosaic the image
                dst = cv2.cvtColor(hduList[0].data, cv2.COLOR_BayerBG2BGR)
                for idx, val in enumerate(dst):
                    hduList[0].data[idx] = val[:,coloridx]
            # Write to temporary file so that we can run solve-field to
            # set WCS data
            hduList.writeto(os.path.join(tmppath, "tmp.fits"), overwrite=True)
            # Now run solve-field to generate final file
            newfits = os.path.join(sciencepath, f)
            rslt = runsolving(hduList[0].header['FOVRA'], hduList[0].header['FOVDEC'],
                os.path.join(tmppath, "tmp.fits"), newfits )
            if rslt == False:
                print("Error solving %s - skipping" % f)
            else:
                tobs = hduList[0].header['MJD-OBS'] * 24 * 60 # MJD in minutes
                # End of accumulator?
                if (len(timeaccumlist) > 0) and ((timeaccumstart + stacktime) < tobs):
                    stackout = os.path.join(stackedpath, "stack-%d.fits" % stackedcnt)
                    tmpout = os.path.join(tmppath, "tmpstack.fits")
                    runstacking(timeaccumra, timeaccumdec, timeaccumlist, stackout)
                    timeaccumlist = []
                    timeaccumstart = 0
                    stackedcnt = stackedcnt + 1
                else:
                    if (timeaccumstart == 0):
                        timeaccumstart = tobs
                        timeaccumra = hduList[0].header['FOVRA']
                        timeaccumdec = hduList[0].header['FOVDEC']
                    timeaccumlist.append(newfits)
                solvedcnt = solvedcnt + 1
            cnt = cnt + 1
    except OSError:
        print("Error: file %s" % f)        
# Final accumulator?
if len(timeaccumlist) > 0:
    stackout = os.path.join(stackedpath, "stack-%d.fits" % stackedcnt)
    runstacking(timeaccumra, timeaccumdec, timeaccumlist, stackout)
    stackedcnt = stackedcnt + 1

print("Processed %d out of %d files and %d stacks into destination '%s'" % (solvedcnt, cnt, stackedcnt, outputdir))

