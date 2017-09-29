# coding=utf-8

import os
import re
import sys

from helpers.helper_config import HelperConfig, FileNotFound
from logger.logger import init_logger

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'


# from clean_mosaic import *
# from combine_mosaics import *
# from combine_surveys import *
# from deconvolve_mosaic import *
# from extraction_hisa import *
# from galprop_skymap import *
# from make_correction import *
# from make_mosaic import *
# from mosaic import *
# from split_mosaic import *


class Survey:
    def __init__(self, survey='MySurvey', species='HI', mosaic='skymap', config_file=False):

        self.logger = init_logger(survey + '_' + mosaic + '_' + species + '_Analysis')
        self.configfilename = 'config/' + survey + '_' + mosaic

        self.helper = HelperConfig(survey_logger=self.logger)

        surveyConfig = self.helper.get_survey_config()
        mosaicConfig = self.helper.get_mosaic_config()
        utilsConfig = self.helper.get_constants_config()
        spectralConfig = self.helper.get_spectral_config()
        spatialConfig = self.helper.get_spatial_config()

        if config_file:
             try:
                 surveyConfigRead, mosaicConfigRead, utilsConfigRead, spectralConfigRead, \
                 spatialConfigRead = self.helper.read_config(self.configfilename)
             except FileNotFound:
                 self.logger.critical("One or more needed files do not exist")
                 return
             try:
                 surveyConfig = self.helper.check_config(surveyConfig, surveyConfigRead)
             except KeyError:
                 return
             try:
                 mosaicConfig = self.helper.check_config(mosaicConfig, mosaicConfigRead)
             except KeyError:
                 return
             try:
                 utilsConfig = self.helper.check_config(utilsConfig, utilsConfigRead)
             except KeyError:
                 return
             try:
                 spectralConfig = self.helper.check_config(spectralConfig, spectralConfigRead)
             except KeyError:
                 return
             try:
                 spatialConfig = self.helper.check_config(spatialConfig, spatialConfigRead)
             except KeyError:
                 return

        self.surveyConf = surveyConfig
        self.mosaicConf = mosaicConfig
        self.utilsConf = utilsConfig
        self.spectralConf = spectralConfig
        self.spatialConf = spatialConfig

        self.flag_existance = False

        self.ret = re.compile('\n')
        self.helper.print(self.surveyConf, 'survey')

    def write_config(self):
        """
        Writes all of the initialization variables to the config file called <surveyname>.cfg.
        """
        self.helper.write_config(self.surveyConf, self.mosaicConf, self.utilsConf, self.spectralConf, self.spatialConf)

    def make_obs(self, type='brightness temperature'):
        """
        Reads the header and gets the data
        """
        try:
            self.obs = Mosaic(self.surveyConf, self.mosaicConf, type)
            self.logger.info(self.ret.subn(', ', str(self.obs))[0])

        except(FileNotFound):
            self.logger.critical("One or more needed files do not exist")
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

        except(FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

    def generate_mosaic(self, species='HI'):
        """
        Generate CGPS-like mosaic.
        Input parameters: species = 'HI'(default),'HISA'(only for CGPS),'HI+HISA'(only for CGPS),
        'CO' (Wco) (only for CGPS)
        Access mosaic's attributes through self.msc
        """
        try:
            self.mosaic
            self.mosaic.species = species
        except AttributeError:
            self.logger.critical("Obs object does not exist. Create it first with the 'makeObs' function.")
            return

        try:
            self.msc = makeMosaic(self.mosaic, self.mosaicConf)
            self.logger.info(self.ret.subn(', ', str(self.msc))[0])
            self.flag_existance = True

        except(FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

    # def loadMosaic(self,species='HI',type='brightness temperature'):
    def load_mosaic(self, species='HI', type='brightness temperature', datatype='original', nmsc=1, totmsc=1,
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
            self.mosaic = Mosaic(self.surveyConf, self.mosaicConf, type, species, datatype, nmsc, totmsc, mypath)
            # self.mosaic = Mosaic(self.surveyConf,self.mosaicConf,type,species,load=True)
            self.logger.info(self.ret.subn(', ', str(self.mosaic))[0])

        except(FileNotFound):
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

        except(FileNotFound):
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
            self.hisa = extractionHISA(self.mosaic, self.spatialConf, self.spectralConf, analysis)
            self.logger.info(self.ret.subn(', ', str(self.hisa))[0])

        except(FileNotFound):
            self.logger.critical("One or more needed files do not exist")
            return

    def get_column_density(self, species='HI'):
        """
        Calculate the column density of a mosaic.
        Keyword parameter: species = 'HI'(default),'HISA'(only for CGPS),'HI+HISA'(only for CGPS)
        Access mosaic's attributes through self.coldens
        """
        try:
            self.mosaic
            self.mosaic.newspec = species
        except AttributeError:
            self.logger.critical(
                "Mosaic object " + species + " does not exist. Create it first with the 'loadMap' function.")
            return
        try:
            self.coldens = makeCorrection(self.mosaic, self.mosaicConf, self.utilsConf)
            self.logger.info(self.ret.subn(', ', str(self.coldens))[0])

        except FileNotFound:
            self.logger.critical("One or more needed files do not exist")
            return

    def combine_mosaics(self, mosaic='skymap', species='HI', type='column density', dim='2D'):

        try:
            self.skyregion = combineMosaics(self.surveyConf, mosaic, species, type, dim)
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
            self.data = plotNvsTs(self.logger, self.obs, self.utilsConf, plot, l, b)
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
        checkForFiles(self.logger, [filename])
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

        except(FileNotFound):
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
            self.mosaic
            self.mosaic.newspec = species
        except AttributeError:
            self.logger.critical(
                "Mosaic object " + species + " does not exist. Load it first with the 'loadMap' function.")
            return
        try:
            self.decon = deconvolveMosaic(self.mosaic, self.mosaicConf, self.utilsConf, rotcurve, scale_data=False)
            self.logger.info(self.ret.subn(', ', str(self.decon))[0])

        except(FileNotFound):
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
        have_mosaic = False
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
                Survey(survey_name, configFile=True)
                print("Analysis start here!!")
                return
            elif opt in ('-i', '--initialize'):
                print("Creating example configuration file called example.cfg")
                mysurvey = Survey(survey_name)
                mysurvey.write_config()
                return

        if not opts:
            raise getopt.GetoptError("Must specify an option, printing help.")

    except getopt.error as e:
        print("Command Line Error: " + e.msg)
        print_cli_help()


if __name__ == '__main__':
    main()
