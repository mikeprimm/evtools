import argparse
from genericpath import isfile
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import os
import logging
from datetime import datetime, timezone
from math import log10, floor
from tenacity import retry, stop_after_delay
import requests

@retry(stop=stop_after_delay(30))
def exofop_getticid(tgtname):
    try:
        url = f"https://exofop.ipac.caltech.edu/tess/gototicid.php?target={tgtname}&json"
        result = requests.get(url)
        rsp = result.json()
        if rsp['status'] == 'OK':
            return rsp['TIC']
        else:
            logging.error(f"EXOFOP error: {rsp['message']}")
            return None
    except Exception as e:
        logging.error(f"EXOFOP error: ${e}")
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
        skycoord = SkyCoord(ra=ra2015*u.deg, 
                          dec=dec2015*u.deg,
                          distance=1*u.mpc, # Far enough to ignore parallax
                          pm_ra_cosdec=pra*u.mas/u.yr,
                          pm_dec=pdec*u.mas/u.yr,
                          obstime=Time('J2015.5'))                
        return skycoord, vmag
    except Exception as e:
        logging.error(f"EXOFOP error: {e}")
        return None, None

