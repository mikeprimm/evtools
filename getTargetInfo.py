import argparse
from genericpath import isfile
import os
import logging
from datetime import datetime, timezone
from math import log10, floor
try:
    from .libs.unistellar import unistellarBestGainAndExp, unstellarExoplanetURL
except ImportError:  # package import
    from libs.unistellar import unistellarBestGainAndExp, unstellarExoplanetURL
try:
    from .libs.exofop import exofop_getcompositeinfo, exofop_getticid, currentRADec
except ImportError: 
    from libs.exofop import exofop_getcompositeinfo, exofop_getticid, currentRADec

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
parser.add_argument("--target", help = "Target Name", required = True)
# Read arguments from command line
try:
    args = parser.parse_args()
except argparse.ArgumentError:
    os.exit(1)

logger.info("args=%s" % args)

tic = exofop_getticid(args.target)

logger.info(f"TIC={tic}")
if tic:
    ra, dec, pra, pdec, vmag = exofop_getcompositeinfo(tic)
    # Adjust for proper motion (milliarcsec to degrees)
    curRA, curDec = currentRADec(ra, dec, pra, pdec)
    logger.info(f"J2000, epoch J2015.5: RA={ra}, Dec={dec}")
    logger.info(f"Current RA={curRA}, Dec={curDec}")
    logger.info(f"V (mag) = {vmag}")
    bestgain, exptime = unistellarBestGainAndExp(vmag)
    logger.info(f"Unistellar: best gain = {bestgain} db, exposure time = {exptime} ms")

url = unstellarExoplanetURL(args.target)
logger.info(f"Unistellar URL (for 1 hour): {url}")
