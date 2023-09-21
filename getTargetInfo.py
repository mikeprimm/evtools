import argparse
from genericpath import isfile
import os
import logging
from astropy.time import Time
from datetime import datetime
from math import log10, floor
try:
    from .libs.unistellar import unistellarBestGainAndExp, unstellarExoplanetURL
except ImportError:  # package import
    from libs.unistellar import unistellarBestGainAndExp, unstellarExoplanetURL
try:
    from .libs.exofop import exofop_getcompositeinfo, exofop_getticid
except ImportError: 
    from libs.exofop import exofop_getcompositeinfo, exofop_getticid

# create logger
logger = logging.getLogger('getTargetInfo')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--target", help = "Target Name", required = True)
parser.add_argument("-d", "--duration", type=float)
# Read arguments from command line
try:
    args = parser.parse_args()
except argparse.ArgumentError:
    os.exit(1)

duration = 1.0
if (args.duration):
    duration = float(args.duration)

logger.info("args=%s" % args)

tic = exofop_getticid(args.target)

logger.info(f"TIC={tic}")
if tic:
    skypos, vmag = exofop_getcompositeinfo(tic)
    # Adjust for proper motion (milliarcsec to degrees)
    curpos = skypos.apply_space_motion(new_obstime = Time(datetime.utcnow(), scale='utc'))
    logger.info(f"J2000, epoch J2015.5: RA/Dec={skypos.to_string('hmsdms')}")
    logger.info(f"Current RA/Dec={curpos.to_string('hmsdms')}")
    logger.info(f"V (mag) = {vmag}")
    bestgain, exptime = unistellarBestGainAndExp(vmag)
    logger.info(f"Unistellar: best gain = {bestgain} db, exposure time = {exptime} ms")

url = unstellarExoplanetURL(args.target, 10)
logger.info(f"Unistellar URL (for 10 seconds): {url}")
url = unstellarExoplanetURL(args.target, duration * 3600)
logger.info(f"Unistellar URL (for {duration} hours): {url}")
