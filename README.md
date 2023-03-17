# Pelagic ocean and climate indicators
This repository includes the code used to create the pelagic ocean and
climate indicators. Each indicator has its own Quarto (.qmd) file which
includes the following:

- Information on accessing the data used for the indicator
- Information on any data wrangling done outside the Quarto file
- Code used to wrangle data and calculate indicator properties
- Code used to generate indicator figures.  Note, though, that most
figures are polished in Adobe Illustrator prior to insertion in the
final report.
- Text that will be copied and pasted into the SAFE report document
  - To the extent possible, this text is formatted to match the SAFE report 
  document.  It's not pretty here, but it simplifies things down the line.
  - Mathematical values related to the indicator (e.g., annual maximum) are
  coded into the text to reduce the risk for failing to update a value.

There's also a `References.qmd` file that includes updated references used 
across all indicators.  `SAFE-Reference-Doc.docx` serves as the style 
template for Quarto.

Note that an additional script is needed for for the Oceanic pH indicator.
This script can be found in the Oceanic_pH folder.

## Questions?  Comments?  Corrections?
Please open an issue or email Phoebe.Woodworth-Jefcoats@noaa.gov

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
