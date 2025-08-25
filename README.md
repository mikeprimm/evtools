# Modified tools for Citizen Science related processing of Exoplanets (original from Mike Primm)

## How to set up
The needed dependencies for these tools are in the provided requirements.txt.  To install these, proceed as normal for your python environment.  Typically this would be done via:

     pip install -r requirements.txt

Python 3.9 or later should work properly, but all code is developed and run on 3.10.

## Available tools

### getTargetInfo.py
This tool can be used to generate deep links for observing an exoplanet.  To run the tool, use the following command:

     python getTargetInfo.py -t "<target-star>" -d <duration-in-hours>

This command will look up the target star (or planet) in EXOFOP, find the coordinates of the star, as well as its V (visual) magnitude, and generate deep links for the requested duration (if not provided, then 1 hour), as well as a link for 10 seconds (for folks who like doing initial darks).

So, to make a deep link for a 3.5 hour observation of TOI-1516b, run:

    python getTargetInfo.py -t "TOI-1516b" -d 3.5

### processExoplanetData.py
This tool is used to both do frame calibration (dark subtraction), color channel extraction, as well as generating substacks for a chosen interval (2 minutes, by default).  Before running the tool, put all darks into one directory, and all science (lights) into another directory.  Then. run the command as follows:

    python processExoplanetData -d <dark-directory-path> -s <science-directory-path> -o <output-directory-path>

By default, this will stack 2 minutes of frames per image (change this with -st \<seconds-to-stack\>), and extract the green channel from the images (control this with -r, -g, -b, -G to select red, green, blue, or grayscale, respectively).   This will also remove any frames that cannot match with the initial frames of each stack (due to clouds or the like), and drop an stacks with less than a minimum of number of matching frames (5 by default, set with -sm \<count\>).

### solveFrames.py
This tool is used to produce frame solutions (sky coordinates) for frames, using a local installation of the astrometry.net command line tools (and associated data) - specifically, the 'solve-field' tool.  This tool isn't strictly necessary for preparation of data for EXOTIC or the like, but does improve results (and can optionally strip any frames where the indicated target is found to be out-of-frame).  To set up the astrometry.net command line tools, and data, see https://g5555.neocities.org/astrometry or other resources describing the process.

Once installed, the (ideally already calbrated and stacked) frames can be used to produce solved frames by running the command:

     python solveFrames.py -i <input-files-path> -o <solved-files-path> -t <target-star-name>

The target star name is optional, but necessary to strip any frames that don't contain the target (due to field drift).  The target lookup is through EXOFOP, so only suitable exoplanet star names will work.

## Typical exoplanet processing steps for using EXOTIC
These are the steps intended for typical use of EXOTIC to do exoplanet analysis, given a raw set of FITS files from the EVScope:

1) Download the files
2) Copy all appropriate dark FITS files into a directory (e.g. darks/)
3) Copy all appropriate science FITS files into another directory (e.g. science/)
4) Run 'python processExoplanetData.py -d <darks-path> -s <science-path> -o <stacks-path>' to build calibrated, 2 minute long substacks
5) (optional) run 'python solveFrames.py -i <stacks-path> -o <final-data-path> -t <exoplanet-name>' to produce solved frames and prune frames not containing the target.
6) Use EXOTIC to process the final frames (just use those frames - no need to deal with debaying or dark subtraction) - how to set up inits.json is outside scope fo these tools.  Remember to use the final set of files - the solved <final-data-path> if step 5 was done, or the <stacks-path> otherwise.
