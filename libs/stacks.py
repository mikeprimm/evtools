from genericpath import isfile
import os
import logging
from astropy.io import fits
import numpy as np

# Build median-based stacked fits HDU from given files
def buildMedianStack(filedir: str, files: list, outfile: str):
    stack = None
    if len(files) > 0:
        stack = fits.HDUList()
        for f in files:
            try:
                dfile = os.path.join(filedir, f)
                # Load file into list of HDU list 
                with fits.open(dfile) as hduList:
                    # Use first one as base
                    if len(stack) == 0:
                        accum = np.zeros((0,) + hduList[0].data.shape)
                        stack.append(hduList[0].copy())
                    accum = np.append(accum, [ hduList[0].data ], axis=0)
            except OSError:
                logging.error("Error: file %s" % f)        
        # Now compute median for each pixel
        accum = np.median(accum, axis=0)
        stack[0].data = accum.astype(np.uint16)
        # And write output, if provided
        if outfile:
            stack.writeto(os.path.join(filedir, outfile), overwrite=True)
    return stack

# Build normalzed stacked flat fits HDU from given files
def buildMasterFlatStack(filedir: str, files: list, outfile: str, darkflat: fits.HDUList):
    normflataccum = None
    if len(files) > 0:
        mflat = fits.HDUList()
        for f in files:
            try:
                dfile = os.path.join(filedir, f)
                # Load file into list of HDU list 
                with fits.open(dfile) as hduList:
                    # Use first one as base
                    if len(mflat) == 0:
                        accum = np.zeros((0,) + hduList[0].data.shape)
                        mflat.append(hduList[0].copy())
                    accum = np.append(accum, [ hduList[0].data ], axis=0)
            except OSError:
                logging.error("Error: file %s" % f)        
        # Now compute median for each pixel
        accum = np.median(accum, axis=0)
        # If we have dark flat, subtract it
        if len(darkflat) > 0:
            # Clamp the data with the dark from below, so we can subtract without rollover
            np.maximum(accum, darkflat[0].data, out=accum)
            # And subtract the dark
            np.subtract(accum, darkflat[0].data, out=accum)
        # And write output dark flat
        mflat[0].data = accum.astype(np.uint16)
        if outfile:
            mflat.writeto(os.path.join(filedir, outfile), overwrite=True)
        # And normalize the flat
        accum = accum.astype(np.float32)
        medi = np.median(accum)
        normflataccum = accum / medi
        # Handle any zero pixels (avoid divide by zero)
        normflataccum[normflataccum == 0] = 1
        #logger.info(f"normflataccum={normflataccum}")
    return normflataccum
