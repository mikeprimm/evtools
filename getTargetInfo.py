import argparse
from genericpath import isfile
import os
import logging
import warnings
from datetime import datetime, timezone
from math import log10, floor
from tenacity import retry, stop_after_delay
import requests

# create logger
logger = logging.getLogger('getTargetInfo')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

@retry(stop=stop_after_delay(30))
def exofop_getticid(tgtname):
    try:
        url = f"https://exofop.ipac.caltech.edu/tess/gototicid.php?target={tgtname}&json"
        result = requests.get(url)
        rsp = result.json()
        if rsp['status'] == 'OK':
            return rsp['TIC']
        else:
            logger.error(f"EXOFOP error: {rsp['message']}")
            return None
    except Exception as e:
        logger.error(f"EXOFOP error: ${e}")
        return None

@retry(stop=stop_after_delay(30))
def exofop_getcompositeinfo(tic):
    try:
        url = f"https://exofop.ipac.caltech.edu/tess/download_target.php?id={tic}"
        result = requests.get(url)
        rsp = result.text.splitlines()
        ra2015 = None
        dec2015 = None
        pra = None
        pdec = None
        vmag = None
        for line in rsp:
            if line.startswith("RA (J2015.5)"):
                sline = line[12:].strip().split(" ")
                ra2015 = float(sline[2])
            elif line.startswith("Dec (J2015.5)"):
                sline = line[13:].strip().split(" ")
                dec2015 = float(sline[2])
            elif line.startswith("Proper Motion RA (mas/yr)"):
                sline = line[25:].strip().split(" ")
                pra = float(sline[0])
            elif line.startswith("Proper Motion Dec (mas/yr)"):
                sline = line[26:].strip().split(" ")
                pdec = float(sline[0])
            elif line.startswith("V     "):
                sline = line[6:].strip().split(" ")
                vmag = float(sline[0])
        return ra2015, dec2015, pra, pdec, vmag
    except Exception as e:
        logger.error(f"EXOFOP error: {e}")
        return None, None, None, None, None

def currentRADec(ra, dec, pra, pdec):
    # Compute position given current time vs J2015.5
    now = datetime.now(timezone.utc)
    j2015_5 = datetime(2015,7,1,0,0,0,0, timezone.utc)
    elapsed_yrs = (now - j2015_5).total_seconds() / 31557600
    # Adjust for proper motion (milliarcsec to degrees)
    curRA = ra + (elapsed_yrs * pra / 3600000.0)
    curDec = dec + (elapsed_yrs * pdec / 3600000.0)

    return curRA, curDec

# Constants from unistellar spreadsheet
baselineExp = 3200.0
baselineGain = 25.0
baselinePeakPixelADU = 3000.0
baselineVmag = 11.7
deltaGain = 5.0
fluxChangeFactor = 1.122 ** deltaGain

def unistellarFluxFromBaseFactor(vmag, exptime):
    return (10 ** ((vmag - baselineVmag)/-2.5)) * (float(exptime)/baselineExp)
def unistellarMaxGain(vmag, exptime):
    return baselineGain - (log10(unistellarFluxFromBaseFactor(vmag, exptime))/log10(1.122))
def unistellarBestGain(vmag, exptime):
    return floor(unistellarMaxGain(vmag, exptime)) - 1
def unistellarBestGainAndExp(vmag):
    exptime = 3970
    bestgain = None
    while exptime >= 1000:
        bestgain = unistellarBestGain(vmag, exptime)
        if (bestgain >= 1):
            return bestgain, exptime
        # Go down by hundreds
        exptime = floor((exptime - 100) / 100) * 100        
    return bestgain, exptime
    
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
