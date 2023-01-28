import argparse
from genericpath import isfile
import os
import pathlib
from astropy.io import fits
import numpy as np
from skimage.color import rgb2gray
from skimage import io
from colour_demosaicing import demosaicing_CFA_Bayer_bilinear

# Initialize parser
parser = argparse.ArgumentParser()
# Add input argument
parser.add_argument("fitsfile", nargs='+', help = "Input files or directories") 
# Adding output argument
parser.add_argument("-o", "--output", help = "Output file or directory") 

# Read arguments from command line
try:
    args = parser.parse_args()
except argparse.ArgumentError:
    os.exit(1)
outputdir='.'
if args.output: 
    outputdir = args.output
filelist=[]
pathlib.Path(outputdir).mkdir(parents=True, exist_ok=True)

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

cnt = 0

for f in filelist:
    try:
        # Load file into list of HDU list 
        with fits.open(f) as hduList:
            basename = os.path.splitext(os.path.basename(f))[0]
            # Demosaic the image
            new_image_data = demosaicing_CFA_Bayer_bilinear(hduList[0].data, "RGGB")
            new_image_data = np.flip(new_image_data, 0)  # Flip Y axis for PNG
            # Normalize to uint8 range (TODO: offer different ranging options vs sqrt)
            dat = np.sqrt((new_image_data - np.min(new_image_data)) / (np.max(new_image_data) - np.min(new_image_data))) * 256
            dat = dat.astype(np.uint8)
            io.imsave(os.path.join(outputdir, "%s.png" % basename), dat, check_contrast=False)
        cnt = cnt + 1
    except OSError:
        print("Error: file %s" % f)        

print("Processed %d files into destination '%s'" % (cnt, outputdir))


