# PySurvey

## Overview

**PySurvey** is the analysis package for studying radio survey data released in FITS format. The analysis is mainly focused on the calulation of column density distributions, and can be performed on mosaics of different surveys (i.e., CGPS, SGPS, VGPS, LAB, Dame) and of two different species, namely atomic hydrogen (HI) and carbon monoxide molecular (CO), a tracer of molecular hydrogen.

To speed up the computation time of the analysis, each mosaic is divided into smaller sub-mosaics. These sub-mosaics can then be processed in parallel if more than one cpu is available. The global variable *glob_ncpu* in **SurveyUtils.py** allows to set the number of cpus that determines the number of sub-mosaics being processed in parrallel. At the end of the analysis all sub-mosaics are brought back together.
