# Pelagic ocean and climate indicators
This repository includes the code used to create the Pelagic Ocean and
Climate Indicators chapter of the SAFE report. In 2026 the code in this repo was
put into a quarto book for easier rendering and to eliminate the copy-paste into
a word document. The Quarto book is rendered to a word document for now. 

The files needed to render the quarto book:
- `_quarto.yaml` is the main configuration file for the quarto book. 
- `_common.R` contains the packages and variables that are used in all 
indicator chapters. *This is where the report year is set and needs to be
updated each year.*
- `index.qmd` contains the introductory text in the chapter, and calls all the 
subsequent indicator chapters. 
- All the chapter indicator quarto (.qmd) files. The most current indicators are:
  - `01-Atmospheric-Carbon-Dioxide.qmd` 
  - `02-Oceanic-pH.qmd`
  - `03-El-Nino-Southern-Oscillation.qmd`
  - `04-Pacific-Decadal-Oscillation.qmd`
  - `05-North-Pacific-Gyre-Oscillation.qmd`
  - `06-Tropical-Cyclones.qmd`
  - `07-Sea-Surface-Temperature.qmd`
  - `08-Temperature-at-Depth.qmd`
  - `09-Ocean-Color.qmd`
  - `10-North_Pacific_Subtropical_Front_(STF)_and_Transition_Zone_Chlorophyll_Front_(TZCF).qmd`
  - `11-Median-Phytoplankton-Size.qmd`
  - `12-Fish-Community-Size-Structure.qmd` Yet to be updated to fit in the Quarto
  book. Will not render as is. 
  - `13-Bigeye_Weight_per_Unit_Effort.qmd` Yet to be updated to fit in the Quarto
  book. Will not render as is. 
  - `14-Bigeye_recruitment_index.qmd` Yet to be updated to fit in the Quarto
  book. Will not render as is. 
  - `15-Bigeye-Forecast.qmd` Yet to be updated to fit in the Quarto book. 
  Will not render as is. 
- `References.bib` contains the BibTex formatted references used across all 
indicators.
- `apa.csl` sets the APA formatting style for the references.
- `SAFE-Reference-Doc-ForBook.docx` serves as the style template for Quarto. 

Each indicator Quarto (.qmd) file includes the following:

- Information on accessing the data used for the indicator
- Information on any data wrangling done outside the Quarto file
- Code used to wrangle data and calculate indicator properties
- Code used to generate indicator text and figures. 
  - There are some checks coded in when indicators change and the text need to
  be updated. These will render as code blocks and need to be commented out
  once the text has been updated. 

Note that an additional script is needed for for the Oceanic pH, Sea Surface
Temperature, Temperature at Depth, Ocean Color, and Median Phytoplankton
Size indicators. The scripts can be found in the eponymous folders.

## Questions? Comments? Corrections?
Please open an issue or email Johanna.Wren@noaa.gov

---

### Disclaimer
This repository is a scientific product and is not official communication 
of the National Oceanic and Atmospheric Administration, or the United 
States Department of Commerce. All NOAA GitHub project code is provided on 
an ‘as is’ basis and the user assumes responsibility for its use. Any 
claims against the Department of Commerce or Department of Commerce bureaus 
stemming from the use of this GitHub project will be governed by all 
applicable Federal law. Any reference to specific commercial products, 
processes, or services by service mark, trademark, manufacturer, or otherwise, 
does not constitute or imply their endorsement, recommendation or favoring by 
the Department of Commerce. The Department of Commerce seal and logo, or the 
seal and logo of a DOC bureau, shall not be used in any manner to imply 
endorsement of any commercial product or activity by DOC or the United 
States Government.

### License
This repository uses the GNU General Public License v3.0 (GPL-3).
Additionally, Software code created by U.S. Government employees 
is not subject to copyright in the United States (17 U.S.C. §105). 
The United States/Department of Commerce reserves all rights to 
seek and obtain copyright protection in countries other than the 
United States for Software authored in its entirety by the Department 
of Commerce. To this end, the Department of Commerce hereby grants 
to Recipient a royalty-free, nonexclusive license to use, copy, and 
create derivative works of the Software outside of the United States.
See LICENSE for further details.
