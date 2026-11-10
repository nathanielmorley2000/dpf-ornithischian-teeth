# Overview

This data and code repository is supplemental to the manuscript:

Morley NED, Wyenberg-Henzler T, Hudgins MN, Sullivan C. 2026. Dental morphology and morphometrics of non-hadrosaurid ornithischian dinosaurs from the Dinosaur Park Formation of Alberta, Canada.

All data and code required to replicate the morphometric analysis are contained within this repository. Scan data for the orodromine tooth is available upon reasonable request from N.E.D. Morley.


# Software Requirements

All .R scripts were created in `R v. 4.5.2`. The following packages are required to execute the code:

- `dplyr`
- `ggord`
- `ggplot2`
- `MASS`
- `openxlsx`
- `pcaMethods`
- `plyr`
- `purrr`
- `readxl`
- `svglite`
- `tibble`
- `tidyr`
- `vegan`

The full references for these packages can be found in Supplementary Table S1 of the corresponding manuscript.


# Repository Guide

### Repository Root

The root of the repository contains the files `dpf-ornithischian-teeth.Rproj`, `.gitignore`, and `README.md`. These files help with repository organization and documentation, and they do not contribute to the results in any way.


### ceratopsidsIncluded

This file contains the data (`dpfornithishcians_ceratopsidsIncluded.xlsx`) and code (`dpfornithischians_ceratopsidsIncluded.R`) required to replicate the morphometric analysis with ceratopsids **included**. Upon execution, the script will create a new subdirectory (`ceratopsidsIncluded/Results`) that will then store the loadings plot and the ordination for the coarse and fine taxonomy CVAs as .svg files, the text output for the CVAs as a .xlsx file, and the text output for the PERMANOVA as a .xlsx file.


### noCeratopsids

This file contains the data (`dpfornithishcians_noCeratopsids.xlsx`) and code (`dpfornithischians_noCeratopsids.R`) required to replicate the morphometric analysis with ceratopsids **excluded**. Upon execution, the script will create a new subdirectory (`noCeratopsids/Results`) that will then store the loadings plot and the ordination for the coarse and fine taxonomy CVAs as .svg files, the text output for the CVAs as a .xlsx file, and the text output for the PERMANOVA as a .xlsx file.


# Contact

If there are any issues with this code, please contact Nathaniel E.D. Morley via email at nmorley@ualberta.ca.

