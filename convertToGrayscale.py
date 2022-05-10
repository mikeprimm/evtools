import argparse
from genericpath import isfile
import os
import pathlib
from astropy.io import fits
import cv2

# Initialize parser
parser = argparse.ArgumentParser()
# Add input argument
parser.add_argument("fitsfile", nargs='+', help = "Input files or directories") 
# Adding output argument
parser.add_argument("-o", "--output", help = "Output file or directory") 
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
filelist=[]
print("output=%s" % outputdir)
if args.fitsfile:
    print("fitsfile=%s" % args.fitsfile)

for f in args.fitsfile:
    if os.path.isdir(f):
        for (dirpath, dirnames, filenames) in os.walk(f):
            for fn in filenames:
                if fn.lower().endswith(".fits"):
                    filelist.append(os.path.join(dirpath,fn))
    elif os.path.isfile(f):
        if f.lower().endswith(".fits"):
            filelist.append(f)
# Make output directory, if needed
pathlib.Path(outputdir).mkdir(parents=True, exist_ok=True)

print("Found %d FITS files" % len(filelist))

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


cnt = 0

for f in filelist:
    try:
        # Load file into list of HDU list 
        with fits.open(f) as hduList:
            basename = os.path.basename(f)
            # DO PROCESSING....
            if togray:
                dst = cv2.cvtColor(hduList[0].data, cv2.COLOR_BayerGB2GRAY)
                for idx, val in enumerate(dst):
                    hduList[0].data[idx] = val
            else:
                # Demosaic the image
                dst = cv2.cvtColor(hduList[0].data, cv2.COLOR_BayerGB2BGR)
                for idx, val in enumerate(dst):
                    hduList[0].data[idx] = val[:,coloridx]
            hduList.writeto(os.path.join(outputdir, basename), overwrite=True)
        cnt = cnt + 1
    except OSError:
        print("Error: file %s" % f)        

print("Processed %d files into destination '%s'" % (cnt, outputdir))


