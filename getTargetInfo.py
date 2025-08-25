import argparse
from genericpath import isfile
import os
import logging
from astropy.time import Time
from datetime import datetime, UTC
from math import log10, floor
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
try:
    from .libs.unistellar import unistellarBestGainAndExp, unstellarExoplanetURL
except ImportError:  # package import
    from libs.unistellar import unistellarBestGainAndExp, unstellarExoplanetURL
try:
    from .libs.exofop import exofop_getcompositeinfo, exofop_getticid, exofop_getparameters
except ImportError: 
    from libs.exofop import exofop_getcompositeinfo, exofop_getticid, exofop_getparameters

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
parser.add_argument("-t", "--target", help = "Target Name")
parser.add_argument("-d", "--duration", type=float)
parser.add_argument("--ra", help = "Target RA")
parser.add_argument("--dec", help = "Target Dec")
parser.add_argument("--mag", help = "Magnitude (V)")

# Read arguments from command line
try:
    args = parser.parse_args()
except argparse.ArgumentError:
    os.exit(1)

duration = 1.0
if (args.duration):
    duration = float(args.duration)
logger.info("args=%s" % args)

if args.target:
    tic = exofop_getticid(args.target)
    logger.info(f"TIC={tic}")
    if tic:
        skypos, vmag = exofop_getcompositeinfo(tic)
        # Adjust for proper motion (milliarcsec to degrees)
        curpos = skypos.apply_space_motion(new_obstime = Time(datetime.now(UTC), scale='utc'))
    else:
        logger.error("Target not found")
        exit(1)
    exofop_getparameters(tic)
elif args.dec and args.ra and args.mag:
    skypos = SkyCoord(args.ra, args.dec, frame='icrs', unit=(u.hourangle, u.deg))
    vmag = float(args.mag)
    curpos = skypos
else:
    logger.error("Either --target or --ra, --dec and --mag required")
    exit(1)

logger.info(f"J2000, epoch J2015.5: RA/Dec={skypos.to_string('hmsdms')}")
logger.info(f"Current RA/Dec={curpos.to_string('hmsdms')}")
logger.info(f"V (mag) = {vmag}")
bestgain, exptime = unistellarBestGainAndExp(vmag)
logger.info(f"Unistellar: best gain = {bestgain} db, exposure time = {exptime} ms")
if bestgain is None:
    exit(1)
elif bestgain == 35:
    logger.warning("Gain is at or exceeding recommended gain ceiling (35dB), so it has been reset to 35dB to prevent saturation. If this target is too faint (>15 Vmag), it may not be visible or require very long stack times")
url = f"unistellar://science/transit?ra={curpos.ra.deg:.5f}&dec={curpos.dec.deg:.5f}&c=3970&et={exptime}&g={bestgain}&d={int(10)}"
logger.info(f"Unistellar URL (for 10 seconds): {url}")
url = f"unistellar://science/transit?ra={curpos.ra.deg:.5f}&dec={curpos.dec.deg:.5f}&c=3970&et={exptime}&g={bestgain}&d={int(duration * 3600)}"
logger.info(f"Unistellar URL (for {duration} hours): {url}")
