# coding=utf-8
import os
from configparser import RawConfigParser


class FileNotFound:
    pass


class CommandNotFound:
    pass


class HelperConfig:
    def __init__(self, survey_logger=None):
        self.logger = survey_logger
        self.constants_dict = {}
        self.survey_dict = {}
        self.mosaic_dict = {}
        self.spectral_dict = {}
        self.spatial_dict = {}

    def set_constants_config(self):
        self.constants_dict = {
            'tcmb': 2.7,            # Cosmic Microwave Background temperature (K)
            'tspin': 125.,          # Excitation or Spin temperature (K) - 125 standard, 150 Fermi
            'xfactor': 1.9e20,      # CO Factor - Strong & Mattox (1996): X=NH2/Wco (K-1 cm-2 km-1 s)
            'c': 1.823e18,          # Costant (cm-2)
            'pc2cm': 3.08567758e18, # Conversion factor from pc to cm (cm)
            'poverk': 4000.,
            'p': 1.0,               # Fraction of HI emission originating behind the HISA cloud
            'fn': 1.0               # Fraction of particle density contributed by the HISA gas, fn = n_hisa/n_tot
        }

    def set_survey_config(self):
        self.survey_dict = {
            'survey': 'MySurvey',
            'species': 'HI'
        }

    def set_mosaic_config(self):
        self.mosaic_dict = {
            'mosaic': 'skymap',
            'lon': None,
            'lat': None,
            'z1': None,
            'z2': None,
            'side': None
        }

    def set_spectral_config(self):
        self.spectral_dict = {
             'n_spectral': 7,           # size of box for spectral smoothing
             'max_loops': 1000,         # maximum number of CLEAN iterations
             'residual_frac': 0.03,     # residual fraction of smoothed for CLEAN cutoff
             'clip_spectral': 0.8,      # fraction of r_max for adding to correction
             'gain_spectral': 0.25,     # fraction of residual height added to correction
             'fwhm_spectral': 8,        # FWHM of the Gaussian for CLEAN loop, in km/s
             'hisa_f_spectral': -2.0,   # residual amplitude factor for potential HISA in spectral search
             'temp_critical': 30.,      # brightness temperature threshold
             'fit_narrow': 2.0,         # FWHM of Gaussian for narrow fit (km/s)
             'fit_broad': 4.0,          # FWHM of Gaussian for broad fit (km/s)
             'fit_qual': 2.0,           # Gaussian fit reliability cutoff
             'dip_cut': 0.6,            # Cutoff for min morphological "dip" for spectral HISA
             'fwhm_spatial_hisa': 5,    # FWHM of Gaussian for spatial smoothing of HISA, in units of pixels
             'min_hisa': 2.0            # cutoff for min HISA amplitude after spatial smoothing
        }

    def set_spatial_config(self):
        self.spatial_dict = {
            'n_spatial': 15,        # size of box for spatial smoothing
            'max_loops': 1000,      # max number of loops for CLEAN algorithm
            'high': 10,             # 10th (or Mth) highest peak in residual used as rmax
            'residual_frac': 0.03,  # fraction of max smoothed height for CLEAN loop cutoff
            'clip_spatial': 0.5,    # fraction of r_max for adding to correction
            'gain_spatial': 0.25,   # fraction of residual added to correction
            'fwhm_spatial': 20,     # FWHM of gaussian for CLEAN, in arcmin
            'noise_resolve': 20,    # angular resolution for calculation of sigma_obs, in minutes
            'hisa_f_spatial': -2.,  # residual amplitude factor for potential HISA in spatial search
            'temp_critical': 30.,   # min unabsorbed temperature for candidate HISA
            'amp_min_first': 4.,    # cutoff amplitude for first HISA filter (pre-smoothing)
            'fwhm_spatial_hisa': 5, # FWHM of Gaussian for spatial smoothing of HISA, in units of pixels
            'min_hisa': 2.          # cutoff for min HISA amplitude after spatial smoothing
        }

    def get_constants_config(self):
        self.set_constants_config()
        return self.constants_dict

    def get_survey_config(self):
        self.set_survey_config()
        return self.survey_dict

    def get_mosaic_config(self):
        self.set_mosaic_config()
        return self.mosaic_dict

    def get_spectral_config(self):
        self.set_spectral_config()
        return self.spectral_dict

    def get_spatial_config(self):
        self.set_spatial_config()
        return self.spatial_dict

    def write_config(self, survey_dict=None, mosaic_dict=None, constants_dict=None, spectral_dict=None, spatial_dict=None):
        """
        Writes all of the needed information to the config file called <mosaicname>.cfg
        """
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

        with open(configfilename + '.cfg', 'w') as configfile:
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
            self.check_for_files([configfilename + ".cfg"])
            self.logger.info('Reading from config file (' + configfilename + '.cfg)')
            config = RawConfigParser()
            config.read(configfilename + '.cfg')

            if config.has_section('survey'):
                survey_dict = dict(config.items('survey'))

            if config.has_section('mosaic'):
                mosaic_dict = dict(config.items('mosaic'))

            if config.has_section('utils'):
                utils_dict = dict(config.items('utils'))

            if config.has_section('spectralSearch'):
                spectral_dict = dict(config.items('spectralSearch'))

            if config.has_section('spatialSearch'):
                spatial_dict = dict(config.items('spatialSearch'))

            return survey_dict, mosaic_dict, utils_dict, spectral_dict, spatial_dict

        except FileNotFound:
            raise FileNotFound

    def check_config(self, referenceDictionary, testDictionary):
        """
        Checks a dictionary against a reference to make sure that all of the parameters are there.
        If all is good, it'll return the checked dictionary. If not, it'll return the reference dictionary
        and raise an exception.
        """
        try:
            for key in referenceDictionary:
                item = testDictionary[key]
            return testDictionary

        except KeyError as inst:
            self.logger.critical("Cannont find " + inst.args[0] + " in the config file.")
            raise KeyError
            return referenceDictionary

    def print(self, configfile, label):
        """
        Prints out information about the various objects to the terminal and to the log file.
        """
        self.logger.info('Reading %s variables...' % label)
        i = 0
        for key, value in configfile.items():
            i += 1
            logString = "%s) %s = %s" % (i, key, value)
            self.logger.info(logString)

    def check_for_files(self, file_list, existence=False):
        """
        Checks for the existence of needed files in the list.
        """
        for filename in file_list:
            if not os.path.exists(filename) and not existence:
                self.logger.critical(filename + " doesn't exist.")
                raise FileNotFound
            elif os.path.exists(filename) and existence:
                self.logger.critical(filename + " already exists.")
                raise FileNotFound

    def checkForCommand(self, command_list):
        """
        Checks for the existence of a certain command.
        """
        for command in command_list:
            cmd = "which -s " + command + " > " + os.devnull + " 2>&1"
            retcode = os.system(cmd)

            if retcode:
                self.logger.critical("unix command " + command + " not found.")
                raise CommandNotFound
