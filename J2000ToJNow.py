import argparse
from datetime import datetime, UTC
from astropy.coordinates import FK4,FK5, SkyCoord
from astropy.time import Time
import astropy.units as u

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("--ra", help = "Target RA")
parser.add_argument("--dec", help = "Target Dec")

# Read arguments from command line
try:
    args = parser.parse_args()
except argparse.ArgumentError:
    os.exit(1)

if (args.ra is None) or (args.dec is None):
   print("Missing --ra or --dec")
   exit(1)
 
# https://astro.swarthmore.edu/transits
# Website reports exoplanet locations in J2000
# Convert them to JNow for use in ASIAir lookup.
 
#_exoCoords = "12:33:03.89 44:54:55.32"  # Copy in RA and Dec from swarthmore with a blank space as a separator.
_exoCoords = "%s %s" % (args.ra, args.dec);


 
# Parse string
_J2000coords = [_exoCoords, _exoCoords.replace(":"," ")]
 
# The logic (Conversion)
_fk5 = FK5(equinox=Time(Time(datetime.now(UTC), scale='utc').jd, format="jd", scale="utc"))
_jNow = SkyCoord(_J2000coords, frame=FK4, unit=(u.hourangle, u.deg)).transform_to(_fk5)
 
# Output
print("\n\nTime now: %s" % datetime.now())
print("J2000 (RA,Dec): %s " % _exoCoords)
result = _jNow[0].to_string('hmsdms').replace("h", "h ")
result = result.replace("m", "m ")
result = result.replace("d", "d ")
print("JNow (RA,Dec): %s \n" % (result))

