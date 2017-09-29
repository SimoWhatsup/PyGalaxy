SURVEY_NAME='MySurvey'
SURVEY_SPECIES='HI'

MOSAIC_NAME='skymap'
MOSAIC_LON='INDEF'
MOSAIC_LAT='INDEF'
MOSAIC_Z1='INDEF'
MOSAIC_Z2='INDEF'
MOSAIC_SIDE='INDEF'

CONSTANT_TCMB=2.7               # Cosmic Microwave Background temperature (K)
CONSTANT_TSPIN=125.             # Excitation or Spin temperature (K) - 125 standard, 150 Fermi
CONSTANT_XFACTOR=1.9e20         # CO Factor - Strong & Mattox (1996): X=NH2/Wco (K-1 cm-2 km-1 s)
CONSTANT_C=1.823e18             # Costant (cm-2)
CONSTANT_PC2CM=3.08567758e18    # Conversion factor from pc to cm (cm)
CONSTANT_POVERK=4000.
CONSTANT_P=1.0                  # Fraction of HI emission originating behind the HISA cloud
CONSTANT_FN=1.0                 # Fraction of particle density contributed by the HISA gas, fn = n_hisa/n_tot

SPECTRAL_NSPECTRAL=7            # size of box for spectral smoothing
SPECTRAL_MAX_LOOPS=1000         # maximum number of CLEAN iterations
SPECTRAL_RESIDUAL_FRAC=0.03     # residual fraction of smoothed for CLEAN cutoff
SPECTRAL_CLIP=0.8               # fraction of r_max for adding to correction
SPECTRAL_GAIN=0.25              # fraction of residual height added to correction
SPECTRAL_FWHM=8                 # FWHM of the Gaussian for CLEAN loop, in km/s
SPECTRAL_HISA_F=-2.0            # residual amplitude factor for potential HISA in spectral search
SPECTRAL_TEMP_CRITICAL=30.      # brightness temperature threshold
SPECTRAL_FIT_NARROW=2.0         # FWHM of Gaussian for narrow fit (km/s)
SPECTRAL_FIT_BROAD=4.0          # FWHM of Gaussian for broad fit (km/s)
SPECTRAL_FIT_QUAL=2.0           # Gaussian fit reliability cutoff
SPECTRAL_DIP_CUT=0.6            # Cutoff for min morphological "dip" for spectral HISA
SPECTRAL_FWHM_SPATIAL_HISA=5    # FWHM of Gaussian for spatial smoothing of HISA, in units of pixels
SPECTRAL_MIN_HISA=2.0           # cutoff for min HISA amplitude after spatial smoothing

SPATIAL_NSPATIAL=15             # size of box for spatial smoothing
SPATIAL_MAX_LOOPS=1000          # max number of loops for CLEAN algorithm
SPATIAL_HIGH=10                 # 10th (or Mth) highest peak in residual used as rmax
SPATIAL_RESIDUAL_FRAC=0.03      # fraction of max smoothed height for CLEAN loop cutoff
SPATIAL_CLIP=0.5                # fraction of r_max for adding to correction
SPATIAL_GAIN=0.25               # fraction of residual added to correction
SPATIAL_FWHM=20                 # FWHM of gaussian for CLEAN, in arcmin
SPATIAL_NOISE_RESOLVE=20        # angular resolution for calculation of sigma_obs, in minutes
SPATIAL_HISA_F=-2.              # residual amplitude factor for potential HISA in spatial search
SPATIAL_TEMP_CRITICAL=30.       # min unabsorbed temperature for candidate HISA
SPATIAL_AMP_MIN_FIRST=4.        # cutoff amplitude for first HISA filter (pre-smoothing)
SPATIAL_FWHM_HISA=5             # FWHM of Gaussian for spatial smoothing of HISA, in units of pixels
SPATIAL_MIN_HISA=2.             # cutoff for min HISA amplitude after spatial smoothing
