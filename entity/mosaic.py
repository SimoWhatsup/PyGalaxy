#!/usr/bin/python
import sys
import numpy as np
from astropy.io import fits
from config import config
from common.util import get_file
from common.logger import get_logger

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'


class Mosaic:
    def __init__(self, survey_conf, mosaic_conf, mtype, species=None, datatype='original', nmsc=0, totmsc=0, path=None):
        """
        Read the fits file: HI, HISA, CO
        To check memory consumption:
            > @profile
            > def function():
            > ...
            > python -m memory_profiler script.py
        :param survey_conf:
        :param mosaic_conf:
        :param mtype: default = brightness_temperature
        :param species:
        :param datatype:
        :param nmsc:
        :param totmsc:
        :param path:
        """
        self.survey = survey_conf['survey']
        self.species = species if species is not None else survey_conf['species']
        self.mosaic = mosaic_conf['mosaic']

        self.mosaic_slug = self.get_mosaic_slug()

        self.logger = get_logger(self.mosaic_slug + '_mosaic')

        self.type = mtype
        self.datatype = datatype
        self.nmsc = nmsc
        self.totmsc = totmsc
        self.path = path

        self.filename, self.mosaic = get_file(self.mosaic_slug, self.datatype, self.nmsc, self.totmsc, self.path)

        # checkForFiles(self.logger, [self.filename])
        self.logger.info('{}'.format(self.filename))

        # Open the file and set the variables
        f = fits.open(self.filename)
        self.keyword = f[0].header

        bscale_flag = False
        if 'bscale' in self.keyword:
            self.bscale = self.keyword['bscale']
            bscale_flag = True
        if 'bzero' in self.keyword and bscale_flag:
            self.bzero = self.keyword['bzero']

        if not ('CROTA1' or 'CROTA2') in self.keyword:
            self.keyword['CROTA1'] = 0.0
            self.keyword['CROTA2'] = 0.0

        # Build arrays
        try:
            self.x, self.y = self.keyword['CRVAL1'], self.keyword['CRVAL2']
            self.dx, self.dy = self.keyword['CDELT1'], self.keyword['CDELT2']
            self.px, self.py = self.keyword['CRPIX1'], self.keyword['CRPIX2']
            self.nx, self.ny = self.keyword['NAXIS1'], self.keyword['NAXIS2']
            self.xarray = self.x + self.dx * (np.arange(self.nx) + 1. - self.px)
            self.yarray = self.y + self.dy * (np.arange(self.ny) + 1. - self.py)
        except:
            self.logger.critical("Some keyword missing. Coordinate arrays cannot be built.")

        if 'ADC_AREA' in self.keyword:
            self.object = self.keyword['ADC_AREA']
            del self.keyword['ADC_AREA']
            self.keyword['OBJECT'] = self.object
        if 'FREQ0' in self.keyword:
            self.band = self.keyword['FREQ0']
            del self.keyword['FREQ0']
            self.keyword['BAND'] = self.band
        elif 'RESTFREQ' in self.keyword:
            self.band = self.keyword['RESTFREQ']
            del self.keyword['RESTFREQ']
            self.keyword['BAND'] = self.band
        elif 'ADC_BAND' in self.keyword:
            self.band = self.keyword['ADC_BAND']
            del self.keyword['ADC_BAND']
            self.keyword['BAND'] = self.band

        if self.keyword['NAXIS'] > 2:
            if not 'CROTA3' in self.keyword:
                self.keyword['CROTA3'] = 0.0
            # Build array
            try:
                self.z = self.keyword['CRVAL3']
                self.dz = self.keyword['CDELT3']
                self.pz = self.keyword['CRPIX3']
                self.nz = self.keyword['NAXIS3']
                self.zarray = self.z + self.dz * (np.arange(self.nz) + 1. - self.pz)
            except:
                self.logger.critical("Some keyword missing. 3rd axis array cannot be built.")

        self.observation = f[0].data
        # free memory
        del f[0]

        self.observation = self.observation.astype(np.float32)

        self.zmin = 1
        self.zmax = self.nz

        if self.type in [config.TB, config.ITB]:
            self.observation[self.observation < -1e4] = 0.
            self.observation[np.isnan(self.observation)] = 0.
            if self.keyword['NAXIS'] > 3:
                if self.survey == 'CGPS':
                    if 'BAND' in self.keyword:
                        if self.keyword['BAND'] == 'HI':
                            self.zmin = 18
                            self.zmax = 271
                            self.observation[:, :self.zmin, :, :] = 0.
                            self.observation[:, self.zmax:, :, :] = 0.
                        if self.keyword['BAND'] == 'CO':
                            self.zmin = 23
                            self.zmax = 256
                            self.observation[:, :self.zmin, :, :] = 0.
                            self.observation[:, self.zmax:, :, :] = 0.
                if self.survey == 'SGPS':
                    self.zmin = 1
                    self.zmax = 410

        self.mosaic = mosaic_conf['mosaic']
        #if not load:
        #    self._inputs = 'Created ' + self.survey + ' Mosaic object ' + self.species + ' ' + self.type
        #else:
        self._inputs = 'Loaded ' + self.survey + ' Mosaic object ' + self.species + ' ' + self.type

    def get_mosaic_slug(self):
        return '{}.{}.{}'.format(self.survey.lower(),
                                 self.species.lower(),
                                 self.mosaic.lower())

    def __repr__(self):
        return self._inputs

    def _obsDialog(self, filename):
        paramDict = MyOrderedDict()
        if filename is None:
            paramDict['filename'] = Param('file', '*.fits')
        else:
            paramDict['filename'] = Param('file', filename)

        root = SimpleDialog(paramDict, title="Mosaic Elements:")
        root.mainloop()
        output = (paramDict['filename'].value())
        return output

    def state(self, output=sys.stdout):
        close = False
        try:
            output = open(output, 'w')
            close = False
        except:
            pass
        output.write("from Mosaic import *\n")
        output.write("obs = MosaicObs(filename=%s)\n" % _quotefn(self.filename))
        if close:
            output.close()
