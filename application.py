#!/usr/bin/python
from common.logger import setup_logging
from survey import Survey

setup_logging()

if __name__ == '__main__':
    # Script to create a mosaic of unabsorbed HI and
    # a corresponding mosaic of column density

    # Canadian Galactic Plane Survey
    cgps = Survey('CGPS', 'HI', 'MW2', False)
    # Create an object (HI or CO, according to the above line)
    cgps.make_obs()
    # Generate mosaics of 'HI' (unabsorbed) or 'HISA' or 'WCO'
    # cgps.generate_mosaic(species='HI')
    # Load new mosaic
    # cgps.load_mosaic(species='HI', mtype='brightness_temperature')
    # Calculate the column density
    cgps.get_column_density(species='HI')
