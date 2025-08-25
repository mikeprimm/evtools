from math import log10, floor
from astropy.time import Time
from datetime import datetime
try:
    from .exofop import exofop_getcompositeinfo, exofop_getticid
except ImportError: 
    from libs.exofop import exofop_getcompositeinfo, exofop_getticid

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
    g = baselineGain - (log10(unistellarFluxFromBaseFactor(vmag, exptime))/log10(1.122))
    if g >= 35:
        return 36
    else:
        return g
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
def unstellarExoplanetURL(target, duration = 3600):
    tic = exofop_getticid(target)
    if tic is None:
        return None
    skypos, vmag = exofop_getcompositeinfo(tic)
    if skypos.ra is None:
        return None
    curpos = skypos.apply_space_motion(new_obstime=Time(datetime.utcnow(), scale='utc'))
    bestgain, exptime = unistellarBestGainAndExp(vmag)
    if bestgain is None:
        return None
    return f"unistellar://science/transit?ra={curpos.ra.deg:.5f}&dec={curpos.dec.deg:.5f}&c=3970&et={exptime}&g={bestgain}&d={int(duration)}"
