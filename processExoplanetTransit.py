import argparse
from genericpath import isfile
import os
import pathlib
from astropy.io import fits
import cv2
import csv
import numpy as np

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

# Read arguments from command line
try:
    args = parser.parse_args()
except argparse.ArgumentError:
    os.exit(1)
outputdir='.'
if args.output: 
    outputdir = args.output
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
    dark.writeto(os.path.join(darkpath, "master=dark.fits"), overwrite=True)

cnt = 0
for f in lightfiles:
    try:
        # Load file into list of HDU list 
        fname = os.path.join(basedir, f)
        with fits.open(fname) as hduList:
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
            hduList.writeto(os.path.join(sciencepath, f), overwrite=True)
            cnt = cnt + 1
    except OSError:
        print("Error: file %s" % f)        

print("Processed %d files into destination '%s'" % (cnt, outputdir))
