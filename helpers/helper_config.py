# coding=utf-8
import os
from configparser import RawConfigParser
from common.logger import get_logger

SURVEY_CONFIG_DIR = 'survey_config/'
logger = get_logger()

class FileNotFound(BaseException):
    pass


class CommandNotFound(BaseException):
    pass


def check_for_files(file_list, existence=False):
    """
    Checks for the existence of needed files in the list.
    """
    for filename in file_list:
        if not os.path.exists(filename) and not existence:
            logger.error(filename + " doesn't exist.")
            raise FileNotFound
        elif os.path.exists(filename) and existence:
            logger.error(filename + " already exists.")
            raise FileNotFound


def check_for_command(command_list):
    """
    Checks for the existence of a certain command.
    """
    for command in command_list:
        cmd = "which -s " + command + " > " + os.devnull + " 2>&1"
        retcode = os.system(cmd)

        if retcode:
            logger.critical("unix command " + command + " not found.")
            raise CommandNotFound


class HelperConfig:
    def __init__(self, name='MySurvey', species='HI', mosaic='skymap'):
        self.logger = get_logger()

        self.name = name
        self.species = species
        self.mosaic = mosaic

        self.constants_config = self.get_constants_config()

        self.survey_config = self.get_survey_config()
        self.mosaic_config = self.get_mosaic_config()
        self.spectral_config = self.get_spectral_config()
        self.spatial_config = self.get_spatial_config()

    def get_constants_config(self):
        return {
            'tcmb': 2.7,  # Cosmic Microwave Background temperature (K)
            'tspin': 125.,  # Excitation or Spin temperature (K) - 125 standard, 150 Fermi
            'xfactor': 1.9e20,  # CO Factor - Strong & Mattox (1996): X=NH2/Wco (K-1 cm-2 km-1 s)
            'c': 1.823e18,  # Constant (cm-2)
            'pc2cm': 3.08567758e18,  # Conversion factor from pc to cm (cm)
            'poverk': 4000.,
            'p': 1.0,  # Fraction of HI emission originating behind the HISA cloud
            'fn': 1.0  # Fraction of particle density contributed by the HISA gas, fn = n_hisa/n_tot
        }

    def get_survey_config(self):
        return {
            'survey': self.name,
            'species': self.species
        }

    def get_mosaic_config(self):
        return {
            'mosaic': self.mosaic,
            'lon': None,
            'lat': None,
            'z1': None,
            'z2': None,
            'side': None
        }

    def get_spectral_config(self):
        return {
            'n_spectral': 7,  # size of box for spectral smoothing
            'max_loops': 1000,  # maximum number of CLEAN iterations
            'residual_frac': 0.03,  # residual fraction of smoothed for CLEAN cutoff
            'clip_spectral': 0.8,  # fraction of r_max for adding to correction
            'gain_spectral': 0.25,  # fraction of residual height added to correction
            'fwhm_spectral': 8,  # FWHM of the Gaussian for CLEAN loop, in km/s
            'hisa_f_spectral': -2.0,  # residual amplitude factor for potential HISA in spectral search
            'temp_critical': 30.,  # brightness temperature threshold
            'fit_narrow': 2.0,  # FWHM of Gaussian for narrow fit (km/s)
            'fit_broad': 4.0,  # FWHM of Gaussian for broad fit (km/s)
            'fit_qual': 2.0,  # Gaussian fit reliability cutoff
            'dip_cut': 0.6,  # Cutoff for min morphological "dip" for spectral HISA
            'fwhm_spatial_hisa': 5,  # FWHM of Gaussian for spatial smoothing of HISA, in units of pixels
            'min_hisa': 2.0  # cutoff for min HISA amplitude after spatial smoothing
        }

    def get_spatial_config(self):
        return {
            'n_spatial': 15,  # size of box for spatial smoothing
            'max_loops': 1000,  # max number of loops for CLEAN algorithm
            'high': 10,  # 10th (or Mth) highest peak in residual used as rmax
            'residual_frac': 0.03,  # fraction of max smoothed height for CLEAN loop cutoff
            'clip_spatial': 0.5,  # fraction of r_max for adding to correction
            'gain_spatial': 0.25,  # fraction of residual added to correction
            'fwhm_spatial': 20,  # FWHM of gaussian for CLEAN, in arcmin
            'noise_resolve': 20,  # angular resolution for calculation of sigma_obs, in minutes
            'hisa_f_spatial': -2.,  # residual amplitude factor for potential HISA in spatial search
            'temp_critical': 30.,  # min unabsorbed temperature for candidate HISA
            'amp_min_first': 4.,  # cutoff amplitude for first HISA filter (pre-smoothing)
            'fwhm_spatial_hisa': 5,  # FWHM of Gaussian for spatial smoothing of HISA, in units of pixels
            'min_hisa': 2.  # cutoff for min HISA amplitude after spatial smoothing
        }

    def write_config(self, survey_dict=None, mosaic_dict=None, constants_dict=None, spectral_dict=None,
                     spatial_dict=None):
        """
        Writes all of the needed information to the config file called <mosaicname>.cfg
        """
        if not os.path.isdir(SURVEY_CONFIG_DIR):
            os.makedirs(SURVEY_CONFIG_DIR)
        configfilename = survey_dict['survey'] + '_' + mosaic_dict['mosaic']

        config = RawConfigParser()
        config.read(configfilename + '.cfg')
        if not config.has_section('survey'):
            config.add_section('survey')

        for variable, value in survey_dict.items():
            config.set('survey', variable, value)
        self.logger.info('wrote common config to ' + configfilename + '.cfg.')

        if mosaic_dict:
            if config.has_section('mosaic'):
                self.logger.info("mosaic config exists, overwriting...")
            else:
                config.add_section('mosaic')
            for variable, value in mosaic_dict.items():
                config.set('mosaic', variable, value)
            self.logger.info("wrote mosaic config to " + configfilename + ".cfg.")

        if constants_dict:
            if config.has_section('constants'):
                self.logger.info("constants config exists, overwriting...")
            else:
                config.add_section('constants')
            for variable, value in constants_dict.items():
                config.set('constants', variable, value)
            self.logger.info("wrote constants config to " + configfilename + ".cfg.")

        if spectral_dict:
            if config.has_section('spectralSearch'):
                self.logger.info("spectralSearch config exists, overwriting...")
            else:
                config.add_section('spectralSearch')
            for variable, value in spectral_dict.items():
                config.set('spectralSearch', variable, value)
            self.logger.info("wrote spectralSearch config to " + configfilename + ".cfg.")

        if spatial_dict:
            if config.has_section('spatialSearch'):
                self.logger.info("spatialSearch config exists, overwriting...")
            else:
                config.add_section('spatialSearch')
            for variable, value in spatial_dict.items():
                config.set('spatialSearch', variable, value)
            self.logger.info("wrote spatialSearch config to " + configfilename + ".cfg.")

        with open(SURVEY_CONFIG_DIR + configfilename + '.cfg', 'w') as configfile:
            config.write(configfile)

    def read_config(self, configfilename):
        """
        Returns all of the needed information from the config file
        called <name>.cfg.  Also checks to make sure all of the
        config parameters are set based on the configure dictionaries
        given in the configDictionaryList.
        """
        survey_dict = {}
        mosaic_dict = {}
        utils_dict = {}
        spectral_dict = {}
        spatial_dict = {}

        try:
            self.check_for_files([SURVEY_CONFIG_DIR + configfilename + ".cfg"])
            self.logger.info('Reading from config file (' + configfilename + '.cfg)')
            config = RawConfigParser()
            config.read(SURVEY_CONFIG_DIR + configfilename + '.cfg')

            if config.has_section('survey'):
                survey_dict = dict(config.items('survey'))

            if config.has_section('mosaic'):
                mosaic_dict = dict(config.items('mosaic'))

            if config.has_section('constants'):
                utils_dict = dict(config.items('constants'))

            if config.has_section('spectralSearch'):
                spectral_dict = dict(config.items('spectralSearch'))

            if config.has_section('spatialSearch'):
                spatial_dict = dict(config.items('spatialSearch'))

            list_of_dict = {
                'survey': survey_dict,
                'mosaic': mosaic_dict,
                'constants': utils_dict,
                'spectral': spectral_dict,
                'spatial': spatial_dict
            }

            return list_of_dict

        except FileNotFound:
            raise FileNotFound

    def check_config(self, reference_dictionary, test_dictionary):
        """
        Checks a dictionary against a reference to make sure that all of the parameters are there.
        If all is good, it'll return the checked dictionary. If not, it'll return the reference dictionary
        and raise an exception.
        """
        try:
            for key in reference_dictionary:
                item = test_dictionary[key]
            return test_dictionary

        except KeyError as inst:
            self.logger.error("Cannot find " + inst.args[0] + " in the config file.")
            raise KeyError

    def print_config(self, configfile, label):
        """
        Prints out information about the various objects to the terminal and to the log file.
        """
        self.logger.info('Reading %s variables...' % label)
        i = 0
        for key, value in configfile.items():
            i += 1
            logString = "%s) %s = %s" % (i, key, value)
            self.logger.info(logString)
