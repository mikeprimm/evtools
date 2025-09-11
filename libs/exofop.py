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
        url = f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic}&json"
        result = requests.get(url)
        rsp = result.json()
        ra2015 = float(rsp['coordinates']['ra']) 
        dec2015 = float(rsp['coordinates']['dec'])
        pra = float(rsp['coordinates']['pm_ra'])
        pdec = float(rsp['coordinates']['pm_dec'])
        vmag = None
        for mag in rsp['magnitudes']:
           if mag['band'] == "V":
               vmag = float(mag['value'])
        skycoord = SkyCoord(ra=ra2015*u.deg, 
                          dec=dec2015*u.deg,
                          distance=1*u.mpc, # Far enough to ignore parallax
                          pm_ra_cosdec=pra*u.mas/u.yr,
                          pm_dec=pdec*u.mas/u.yr,
                          obstime=Time('J2015.5'))   
        return skycoord, vmag
    except Exception as e:
        print(f"ERROR %s" % e)
        logging.error(f"EXOFOP error: {e}")
        return None, None

@retry(stop=stop_after_delay(30))
def exofop_getparameters(tic):
    try:
        url = f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic}&json"
        result = requests.get(url)
        rsp = result.json()
        print(rsp.keys())

        # print(rsp['coordinates'])
        print(rsp['planet_parameters'])
        # print(rsp['stellar_parameters'])
        return None, None
    except Exception as e:
        logging.error(f"EXOFOP error: {e}")
        return None, None
