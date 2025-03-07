---
id: config
title: Configuration Files
---

 The configuration files describe an individual instrument and contain all the information for the physics-based fundamental parameters model. This information is described in Section 6. Note that the information in this file must be accurate to get reliable calculations or quantification.

Configuration files are generally identical to the spectrum files with the same file extension but with zero number of spectrum data points. Several of the keywords in the ISO standards that these file formats are based on are required in the standard. These keywords are included in the configuration files even though they are not relevant to the configuration but only to the spectrum data. A spectrum can be included in a configuration file but will be ignored.

## Old XSP Configuration Files

These files are an older format that was used for the Mars Borehole XRF instrument. They were used to analyze data for the PIXL breadboard instrument before the PIXL instrument was selected and during the early design process. It was based on an earlier version of the ISO standard in the next section. That older standard did not provide for user-defined keywords so several ad-hoc keywords were created for this format. This format is gradually being replaced by the new MSA configuration file format.

#FORMAT : EMSA/MAS Spectral Data File
#VERSION : 1.0
#TITLE : PIXL best model for current FM configuration (1 hr, side window tube, 32 mm standoff, 2 det)
#DATE : 07/01/2015
#TIME : 0:00:00 AM
#OWNER : Tim
#NPOINTS : 0
#NCOLUMNS : 1
#XUNITS : eV
#YUNITS : COUNTS
#DATATYPE : Y
#XPERCHAN : 10
#OFFSET : 0
#LIVETIME : 3600
#SIGNALTYPE : XRF
anode_z: 45
BEAMKV: 28.00
tube_inc_angle: 30
tube_takeoff_angle: 70
tube_be_window: 0.150
tube_current: 0.02
filter_z: 1
filter_thick: 0
excit_angle: 90.00
emerg_angle: 70.00
solid_angle: 0.00000912
path_type: 2
inc_path_length: 3.0
emerg_path_length: 3.2
window_type: 2
window_thick: 0
detector_type: 2
minimum_energy: 900
optic_type: 3
#SPECTRUM
#ENDOFDATA

## New MSA Configuration Files

These files follow the ISO 22029 standard†, which details the EMSA/MAS standard file format for spectral data exchange. The format is a keyword and value format and was developed for X-ray emission spectra taken with electron microscopes. User-defined keywords are provided in the standard and several have been defined for PIXL (and are also applicable to other XRF instruments). The format has also been adapted to accommodate multiple detectors. Multiple columns are allowed in the standard for spectrum data and this was extended to include multiple entries for energy calibration, live time, and similar keywords.

#FORMAT : EMSA/MAS spectral data file
#VERSION : TC202v2.0 PIXL
#COMMENT : The first two lines must appear exactly as shown, other keywords can be in any order.
#COMMENT : Updated Feb. 7, 2018 to include ##TRIGGERS and ##EVENTS keywords
#TITLE : This is a descriptive title (NB: this keyword must be present even if no title)
#TITLE : This file is based on the PIXL Breadboard configuration of August 2017
#TITLE : The title keyword can be repeated several times if necessary
#COMMENT : Comments can appear anywhere except in the spectrum data.
#DATE : Date in the format DD-MMM-YYYY, for example 30-SEP-2017
#TIME : 12:22 The time of day at which the spectrum was recorded, in 24-hour format
#OWNER : NewBB (PIXL will use this to indicate which instrument took the data.)
#NPOINTS : 0 This should be zero for configuration only files
#NCOLUMNS : 2 This will be ignored for configuration files
#XUNITS : eV
#YUNITS : COUNTS
#DATATYPE : Y (This would be YY for two detectors.)
#XPERCHAN : 10.0, 10.0 eV per channel, will be used if the spectrum file is not calibrated
#OFFSET : 0.00, 0.0 eV of first channel, will be used if the spectrum file is not calibrated
#SIGNALTYPE : XRF
#COMMENT : The above keywords (except the comments) are always required in the ISO 22029:2012(E) format. Generally, the keywords that apply to all spectra will be in the configuration file and only the keywords for information unique to each spectrum will appear in the spectrum files. However, any keywords can appear in the spectrum files and will override the configuration values. The keywords below are instrument configuration information and will normally by in the configuration file.
##ANODE : 45 Atomic number of anode in X-ray tube
#BEAMKV : 28.00 X-ray tube voltage in kilovolts
##TUBEINCANG : 90.00 X-ray tube electron incident angle in degrees
##TUBETAKEOF : -90.0 X-ray tube takeoff angle in degrees (negative for transmission anode)
##TUBEWINDOW : 0.275 X-ray tube Be window thickness in mm
#EMISSION : 20 X-ray tube emission current in microAmps
##FILTERZ : 1 Primary beam filter material - atomic number of metal foil
##FILTERTH : 0 Primary beam filter thickness in microns
##OPTICFILE : 5 File name of optic transmission function (5 for 2017 breadboard, 3 for old PIXL breadboard, 0 for no optic)
##INCSR : 0.0017 Solid angle from source in steradians (can include normalization for optic file)
##INCANGLE : 90.00 Incident angle of primary X-ray beam in degrees (90 is normal incidence)
#ELEVANGLE : 70.00 Elevation angle of detector, in degrees (90 is normal to surface)
#AZIMANGLE : 180.0 Azimuth angle between incident beam plane and detected beam plane
##GEOMETRY : 1.0 Geometric correction factor
#SOLIDANGLE : 1.571 Solid angle collected by the detector in steradians
#COMMENT : This is for one detector, double for two detectors
#EDSDET : SDBEW Type of XRF detector (SDBEW, SIBEW, CDBEW, or GEBEW, for SDD, Si_PIN, CdTe, or Ge)
#TBEWIND : 0.0017 Thickness of Be window on detector in cm
#TACTLYR : 0.050 Thickness of active layer of detector in cm
##DETRES : 129 Detector energy resolution in eV (at 5.9 keV, Mn Ka emission line, Fe55 source)
##ATMOSPHERE : He Atmosphere in X-ray beam path, can be Vac, He, Mars, Earth, Air
##PATHINCLEN : 2.0 Length of incident beam path in cm
##PATHEMGLEN : 3.2 Length of emerging beam path in cm
##WINDOWTYPE : None Type of window between instrument and specimen (None, B4C, Plastic, Zr, Al, Nylon, or Al2O3)
##WINDOWTH : 0.00 Thickness of above window in microns (0 = No window)
#COMMENT : Live time is given in the configuration file for calculations without reading a measured spectrum
#LIVETIME : 1.0, 1.0 seconds (separate by commas if more than one detector
##MINIMUM_EN : 900 Minimum energy in eV for any emission lines to be included in the spectrum
#COMMENT : Optional keywords we may want to consider in the future:
#COMMENT : #CHECKSUM (See definition in standard)
#COMMENT : #BEAMDIAM beam diameter on specimen (nanometers)
#COMMENT : #THICKNESS thickness of specimen (nanometers)

† ISO 22029:2012. Microbeam analysis — EMSA/MAS standard file format for spectral-data exchange. See International Organization for Standardization, https:// www.iso.org/standard/56211.html 