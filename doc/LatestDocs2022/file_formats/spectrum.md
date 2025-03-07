---
id: spectrum
title: Spectrum Files
---

 ## MSA Format

This format is identical to the new MSA configuration files. The number of data points must be correct for the spectrum data. The spectrum data is read until the specified number of points has been read. If an end of file or the end of data marker are found before all of the points are read an error results. Many of the configuration keywords are usually omitted from the spectrum files but all of the required keywords from the relevant ISO standard must be present.

An example of an MSA spectrum file from the PIXL Breadboard instrument is given here with most of the actual spectrum data removed. This format is produced when maps are collected using the LabView software on the PIXL breadboard or the Stony Brook instrument. A few comments have been added that are not written by the LabView software.

#FORMAT : EMSA/MAS spectral data file
#VERSION : TC202v2.0 PIXL
#COMMENT : Updated Feb. 7, 2018 to include ##TRIGGERS and ##EVENTS keywords
#TITLE : Control Program v6 Test 5
#OWNER : JPL BREADBOARD vx
#DATE : 01-29-2018
#TIME : 09:27:54
#NPOINTS : 4096
#NCOLUMNS : 2
#XUNITS : eV
#YUNITS : COUNTS
#DATATYPE : YY
#XPERCHAN : 10.0, 10.0 eV per channel
#OFFSET : 0.0, 0.0 eV of first channel
#SIGNALTYPE : XRF
#COMMENT : Comments like this can be included anywhere except in the spectrum data.
#COMMENT : Fields above will be written to config file and spectrum files; fields below only config file.
#XPOSITION : 0.000
#YPOSITION : 0.000
#ZPOSITION : 0.000
#LIVETIME : 121.0, 121.0
#REALTIME : 121.1, 121.1
#COMMENT : The livetime from XIA pulse processing is not the usual livetime. It must be adjusted
#COMMENT : by the ratio of output to input counts. These quantities are read from the pulse processor and included
#COMMENT : using the following two keywords. PIQUANT will make the correction if it finds these keywords.
##TRIGGERS : 194764, 195575
##EVENTS : 190020, 190882
#SPECTRUM :
0, 0
0, 0
0, 0
0, 0
#ENDOFDATA :
 
## MCA Format

These spectrum files are specific to the software that is provided with the Amptek and Ketek X-ray detectors. Both detector formats can be read and are distinguished by reading the first line of the file. The format is not described here and was discovered empirically from files written by the respective detector software. This format is not written by any of the PIXL instruments unless the software from the detector manufacturer is being used to collect spectra.

## XSP Format

This format is identical to the old XSP configuration files. The number of data points must be correct for the spectrum data. The spectrum data is read until the specified number of points has been read. If an end of file or the end of data marker are found before all the points are read an error results. Some of the configuration keywords can be omitted from the spectrum files but all the required keywords from the old ISO standard must be present. This format is not currently written by any PIXL instruments or software. 