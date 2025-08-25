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

def safe_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

def _get_vmag(data):
    for m in data.get("magnitudes", []):
        band = (m.get("band") or "").strip().lower()
        if band == "v":
            v = safe_float(m.get("value"))
            if v is not None:
                return v
    return None

@retry(stop=stop_after_delay(30))
def exofop_getcompositeinfo(tic):
    try:
        url = f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic}&json"
        result = requests.get(url, timeout=15)
        result.raise_for_status()
        rsp = result.json()
        coords = rsp.get("coordinates", {})
        ra2015 = safe_float(coords.get("ra"))
        dec2015 = safe_float(coords.get("dec"))
        pra = safe_float(coords.get("pm_ra"))
        pdec = safe_float(coords.get("pm_dec"))
        if None in (ra2015, dec2015, pra, pdec):
            raise ValueError("Missing or invalid RA/Dec/PM fields in EXOFOP JSON.")
        dist_pc = None
        for sp in rsp.get("stellar_parameters", []):
            dist_pc = safe_float(sp.get("dist"))
            if dist_pc is not None:
                break
        distance = (dist_pc * u.pc) if dist_pc is not None else (1 * u.mpc)
        vmag = _get_vmag(rsp)
        skycoord = SkyCoord(
            ra=ra2015 * u.deg,
            dec=dec2015 * u.deg,
            distance=distance,
            pm_ra_cosdec=pra * u.mas / u.yr,
            pm_dec=pdec * u.mas / u.yr,
            obstime=Time("J2015.5")
        )
        return skycoord, vmag
    except Exception as e:
        logging.error(f"EXOFOP error for TIC {tic}: {e}")
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
