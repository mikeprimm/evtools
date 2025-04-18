from genericpath import isfile
import os
import logging
from astropy.io import fits
import numpy as np

# Fix raw data if needed (workaround for bad Unistellar export)
def getDataBuffer(hduList: fits.HDUList):
    data = hduList[0].data
    if 'SOFTVER' in hduList[0].header:
        if hduList[0].header['SOFTVER'] == '4.1-6ded42a5':
            np.multiply(data, 2, out=data)
    return data

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
                    data = getDataBuffer(hduList)
                    # Use first one as base
                    if len(stack) == 0:
                        accum = np.zeros((0,) + data.shape)
                        stack.append(hduList[0].copy())
                    accum = np.append(accum, [ data ], axis=0)
            except OSError:
                logging.error("Error: file %s" % f)        
        # Now compute median for each pixel
        accum = np.median(accum, axis=0)
        stack[0].data = accum.astype(np.uint16)
        # Handle broken size for early Odyssey darks
        if stack[0].data.shape[0] == 2190:
            print("Resize broken Odyssey darks")
            newsize = np.zeros((2192, stack[0].data.shape[1]), dtype=np.uint16)
            newsize[1:2191] = stack[0].data[0:2190]
            stack[0].data = newsize
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
                    data = getDataBuffer(hduList)
                    # Use first one as base
                    if len(mflat) == 0:
                        accum = np.zeros((0,) + data.shape)
                        mflat.append(hduList[0].copy())
                    accum = np.append(accum, [ data ], axis=0)
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
        #medi = [ [ np.median(accum[::2,::2]), np.median(accum[::2,1::2])],
        #         [ np.median(accum[1::2,::2]), np.median(accum[1::2,1::2])] ]
        #medi = np.tile(medi, (accum.shape[0]//2,accum.shape[1]//2))
        medi = np.median(accum)
        normflataccum = accum / medi
        # Handle any zero pixels (avoid divide by zero)
        normflataccum[normflataccum == 0] = 1.0
        #logger.info(f"normflataccum={normflataccum}")
    return normflataccum

# Scale single color channel up 4x4 for stacking
#   Each pixel is treated as a 4x4 grid centered
#   For red/blue (one pixel ever other row and column), the corresponding quarters of the intermediate
#   pixels are assumed to match the pixel:
#
#        ---              ------------
#        -R-    becomes   ------------
#        ---              --RRRRRRRR--
#                         --RRRRRRRR--
#                         --RRRRRRRR--
#                         --RRRRRRRR--
#                         --RRRRRRRR--
#                         --RRRRRRRR--
#                         --RRRRRRRR--
#                         --RRRRRRRR--
#                         ------------
#                         ------------
#
#   For green (2 pixels out of each 2x2), the nearest 25% of the intermediate pixels are assumed
#   to match (diamond pattern):
#
#        ---              ------------
#        -G-    becomes   ------------
#        ---              -----G------
#                         ----GGG-----
#                         ----GGGGG---
#                         ---GGGGGGG--
#                         --GGGGGGG---
#                         ---GGGGG----
#                         -----GGG----
#                         ------G-----
#                         ------------
#                         ------------
#
def scaleAndDemosaicImage(data: np.array):
    print(f"data.shape={data.shape}")
    newshape = (data.shape[0]*4, data.shape[1]*4)
    clrshape = (data.shape[0]>>1, data.shape[1]>>1)
    reddata = np.zeros(clrshape, dtype=data.dtype)
    green1data = np.zeros(clrshape, dtype=data.dtype)
    green2data = np.zeros(clrshape, dtype=data.dtype)
    bluedata = np.zeros(clrshape, dtype=data.dtype)
    redxmap = np.minimum((np.arange(newshape[0]) + 2) >> 3, clrshape[0]-1)
    redymap = np.minimum((np.arange(newshape[1]) + 2) >> 3, clrshape[1]-1)
    reddata[:,:] = data[::2,::2]       
    green1data[:,:] = data[::2,1::2]
    green2data[:,:] = data[1::2,::2]
    bluedata[:,:] = data[1::2,1::2]
    red = np.zeros(newshape, dtype=data.dtype)
    green = np.zeros(newshape, dtype=data.dtype)
    blue = np.zeros(newshape, dtype=data.dtype)
    for x in range(len(redxmap)):
        red[x,:] = reddata[redxmap[x],redymap]
        blue[x,:] = bluedata[redxmap[x],redymap]
    print(f"red={red}")
    print(f"blue={blue}")
    return red, green, blue