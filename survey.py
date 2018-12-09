# coding=utf-8
import os
import re
import sys

from helpers.helper_config import HelperConfig, FileNotFound, check_for_files
from common.logger import get_logger

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from entity.mosaic import Mosaic
# from clean_mosaic import *
# from combine_mosaics import *
# from combine_surveys import *
# from deconvolve_mosaic import *
# from extraction_hisa import *
# from galprop_skymap import *
from make_correction import MakeCorrection
from make_mosaic import MakeMosaic
# from split_mosaic import *


class Survey:
    def __init__(self, survey='MySurvey', species='HI', mosaic='skymap', read_config=False):

        self.logger = get_logger(survey + '_' + mosaic + '_' + species + '_Analysis')
        self.configfilename = survey + '_' + mosaic

        self.helper = HelperConfig(name=survey, species=species, mosaic=mosaic)

        survey_config = self.helper.survey_config
        mosaic_config = self.helper.mosaic_config
        constants_config = self.helper.constants_config
        spectral_config = self.helper.spectral_config
        spatial_config = self.helper.spatial_config

        if read_config:
            try:
                config_read_dict = self.helper.read_config(self.configfilename)

                survey_config_read = config_read_dict.get('survey')
                mosaic_config_read = config_read_dict.get('mosaic')
                constants_config_read = config_read_dict.get('constants')
                spectral_config_read = config_read_dict.get('spectral')
                spatial_config_read = config_read_dict.get('spatial')

            except FileNotFound:
                self.logger.error('One or more needed files do not exist')
                return

            survey_config = self.helper.check_config(survey_config, survey_config_read)
            mosaic_config = self.helper.check_config(mosaic_config, mosaic_config_read)
            constants_config = self.helper.check_config(constants_config, constants_config_read)
            spectral_config = self.helper.check_config(spectral_config, spectral_config_read)
            spatial_config = self.helper.check_config(spatial_config, spatial_config_read)

        self.survey_conf = survey_config
        self.mosaic_conf = mosaic_config
        self.utils_conf = constants_config
        self.spectral_conf = spectral_config
        self.spatial_conf = spatial_config

        self.flag_existance = False

        self.ret = re.compile('\n')
        self.helper.print_config(self.survey_conf, 'survey')

    def write_config(self):
        """
        Writes all of the initialization variables to the config file called <surveyname>.cfg.
        """
        self.helper.write_config(self.survey_conf, self.mosaic_conf, self.utils_conf, self.spectral_conf, self.spatial_conf)

    def make_obs(self, mtype='brightness_temperature'):
        """
        Reads the header and gets the data
        """
        try:
            self.mosaic = Mosaic(self.survey_conf, self.mosaic_conf, mtype)
            self.logger.info(self.mosaic)

        except FileNotFound:
            self.logger.error('One or more needed files do not exist')
            return

    def clean_mosaic(self, scale_data=False):
        """
        Access mosaic's attributes through self.clean
        Input parameters:
            species = 'HI'(default), 'CO'
            scale   = True (write scaled data in FITS file), False (do not scale data)
        """
        try:
            self.mosaic
        except AttributeError:
            self.logger.critical("Mosaic object does not exist. Create it first with the 'loadMosaic' function.")
            return

        try:
            self.clean = CleanMosaic(self.mosaic, scale_data)
            self.logger.info(self.ret.subn(', ', str(self.clean))[0])

        except FileNotFound:
            self.logger.critical("One or more needed files do not exist")
            return

    def generate_mosaic(self, species='HI'):
        """
        Generate CGPS-like mosaic.
        Input parameters: species = 'HI'(default), 'HISA'(only for CGPS), 'HI+HISA'(only for CGPS),
        'CO' (Wco) (only for CGPS)
        Access mosaic's attributes through self.msc
        """
        try:
            self.mosaic.species = species
        except AttributeError:
            self.logger.error('Obs object does not exist. Create it first with the \'make_obs\' function.')
            return

        try:
            self.msc = MakeMosaic(self.mosaic, self.mosaic_conf)
            self.logger.info(self.ret.subn(', ', str(self.msc))[0])
            self.flag_existance = True

        except FileNotFound:
            self.logger.error('One or more needed files do not exist')
            return

    def load_mosaic(self, species='HI', mtype='brightness_temperature', datatype='original', nmsc=1, totmsc=1,
                    mypath=None):
        """
        Load a mosaic.
        Input parameters:
            - species  = 'HI'(default),'HISA'(only for CGPS),'HI+HISA'(only for CGPS column density),
                    'CO' (Wco) (only for CGPS)
            - type     = 'brightness temperature'(defualt),'column density'
            - datatype = 'original' (default),'clean' (after applying data-clean methods, see cleanMosaic)
            - load     = True, False (original survey mosaic)
        Access mosaic's attributes through self.mosaic
        """
        try:
            self.mosaic = Mosaic(self.survey_conf, self.mosaic_conf, mtype, species, datatype, nmsc, totmsc, mypath)
            # self.mosaic = Mosaic(self.survey_conf,self.mosaic_conf,type,species,load=True)
            self.logger.info(self.ret.subn(', ', str(self.mosaic))[0])

        except FileNotFound:
            return

    def get_galprop_map(self, resolution=1):
        """
        Access mosaic's attributes through self.galprop
        Input parameters:
            - resolution, in degrees
        """
        try:
            self.mosaic
        except AttributeError:
            self.logger.critical("Mosaic object does not exist. Create it first with the 'loadMosaic' function.")
            return

        try:
            self.galprop = GalpropSkymap(self.mosaic, resolution)
            self.logger.info(self.ret.subn(', ', str(self.galprop))[0])

        except FileNotFound:
            return

    # spectral search needs to be implemented
    def extract_hisa(self, analysis='spatial'):
        """
        Access mosaic's attributes through self.hisa
        """
        try:
            self.mosaic
        except AttributeError:
            self.logger.critical(
                "Mosaic object " + species + " does not exist. Create it first with the 'loadMap' function.")
            return

        try:
            self.hisa = extractionHISA(self.mosaic, self.spatial_conf, self.spectral_conf, analysis)
            self.logger.info(self.ret.subn(', ', str(self.hisa))[0])

        except FileNotFound:
            self.logger.critical("One or more needed files do not exist")
            return

    def get_column_density(self, species='HI'):
        """
        Calculate the column density of a mosaic.
        Keyword parameter: species = 'HI'(default),'HISA'(only for CGPS),'HI+HISA'(only for CGPS)
        Access mosaic's attributes through self.coldens
        """
        try:
            self.mosaic.newspec = species
        except AttributeError:
            self.logger.critical(
                "Mosaic object " + species + " does not exist. Create it first with the 'loadMap' function.")
            return
        try:
            self.coldens = MakeCorrection(self.mosaic, self.mosaic_conf, self.utils_conf)
            self.logger.info(self.ret.subn(', ', str(self.coldens))[0])

        except FileNotFound:
            self.logger.critical("One or more needed files do not exist")
            return

    def combine_mosaics(self, mosaic='skymap', species='HI', type='column density', dim='2D'):
        try:
            self.skyregion = combineMosaics(self.survey_conf, mosaic, species, type, dim)
            self.logger.info(self.ret.subn(', ', str(self.skyregion))[0])

        except FileNotFound:
            self.logger.critical("One or more needed files do not exist")
            return

    def combine_surveys(self, surveylist, mosaiclist, species='HI', res=0.125):

        try:
            self.sky = combineSurveys(surveylist, mosaiclist, species, res)
            self.logger.info(self.ret.subn(', ', str(self.sky))[0])

        except FileNotFound:
            self.logger.critical("One or more needed files do not exist")
            return

    def make_plot(self, plot='NH vs Ts', l=0., b=0.):

        try:
            self.data = plotNvsTs(self.logger, self.obs, self.utils_conf, plot, l, b)
            self.logger.info(self.ret.subn(', ', str(self.data))[0])

        except FileNotFound:
            self.logger.critical("One or more needed files do not exist")
            return

    def delete_mosaic(self):
        """
        Delete the last object loaded with the loadMosaic function.
        Usage:
            survey.loadMosaic(species='some',type='some')
            deleteMosaic()
        """
        try:
            self.mosaic
        except AttributeError:
            self.logger.critical("Cannot delete Mosaic object because does not exist.")
            return

        filename = self.mosaic.mosaicFile
        check_for_files([filename])
        self.logger.info('Removing file: ' + filename)
        os.remove(filename)

    def split_mosaic(self, num=1):
        """
        Access mosaic's attributes through self.test
        """
        try:
            self.mosaic
        except AttributeError:
            self.logger.critical(
                "Mosaic object " + species + " does not exist. Load it first with the 'loadMap' function.")
            return
        try:
            self.split = splitMosaic(self.mosaic, num)
            self.logger.info(self.ret.subn(', ', str(self.split))[0])

        except FileNotFound:
            self.logger.critical("One or more needed files do not exist")

    def deconvolve_mosaic(self, species='HI', rotcurve='Bissantz2003'):
        """
        Deconvolution technique based on the Galactic rotation curve.
        Keyword parameters:
            species
            rotcurve = 'Bissantz2003', 'Clemens1985'
        Access mosaic's attributes through self.decon
        """
        try:
            self.mosaic.newspec = species
        except AttributeError:
            self.logger.critical(
                "Mosaic object " + species + " does not exist. Load it first with the 'loadMap' function.")
            return
        try:
            self.decon = deconvolveMosaic(self.mosaic, self.mosaic_conf, self.utils_conf, rotcurve, scale_data=False)
            self.logger.info(self.ret.subn(', ', str(self.decon))[0])

        except FileNotFound:
            return


def print_cli_help():
    """
    This function prints out the help for the CLI.
    """
    cmd = os.path.basename(sys.argv[0])

    print(
        """
                       - Survey - 
           
           Perform Survey analysis for different line emissions. 
           You can use the command line functions listed below or run this module
           from within python.
           
           %s (-h|--help) ... This help text.
           
           %s (-a|--analyze) (-n |--surveyname=)<surveyname> ...  Perform an analysis
           on <surveyname>.  <surveyname> is the prefix used for this analysis.
           You must already have a configuration file if using the command
           line interface.
           
           %s (-i|--initialize) ... Generate a default config file called
           example.cfg.  Edit this file and rename it <surveyname>.cfg for use
           in the Survey module.
           
           """ % (cmd, cmd, cmd))


def main():
    """
    Command-line interface.  Call this without any options for usage notes.
    """
    import getopt
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hiamxcb:n:', ['help', 'analyze', 'initialize', 'surveyname', ])
        # Loop through first and check for the surveyname
        have_mosaic = True  # False
        survey_name = 'example'
        for opt, val in opts:
            if opt in ('-n', '--surveyname'):
                have_mosaic = True
                survey_name = val
            elif opt in ('-h', '--help'):
                print_cli_help()
                return
            elif opt in ('-a', '--analyze'):
                if not have_mosaic:
                    raise getopt.GetoptError("Must specify surveyname, printing help.")
                Survey(survey=survey_name, read_config=True)
                print("Analysis start here!!")
                return
            elif opt in ('-i', '--initialize'):
                print("Creating example configuration file called example.cfg")
                mysurvey = Survey(survey=survey_name)
                mysurvey.write_config()
                return

        if not opts:
            raise getopt.GetoptError("Must specify an option, printing help.")

    except getopt.error as e:
        print("Command Line Error: " + e.msg)
        print_cli_help()


if __name__ == '__main__':
    main()
