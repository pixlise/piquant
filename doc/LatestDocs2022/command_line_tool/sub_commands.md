---
id: sub_commands
title: Sub-Commands
---

 Each of the actions described in the GUI sections has a corresponding sub-command in the command line tool. The GUI sections have a more complete description of how each sub-command works and what the arguments mean, especially the element fit control list.

Only the first three letters of the sub-command need to be entered except for calculate, which needs four letters. The sub-commands are not case sensitive. The full word can be entered but only the first three or four letters are checked. The sub-command cannot contain any blanks.

Only a brief description, the expected arguments, and any special notes on use of the command line are included here. See the corresponding sections in the GUI actions for more general usage information.

## Energy Calibration

Associates the two largest peaks with elements from the element fit control list and calculates the energy of the first channel and the energy per channel.

### Argument list:

* Spectrum file (read)
* Element fit control list (one or two elements)

## Plot

This sub-command converts a spectrum file from the accepted formats to a comma-separated-values file with energy in the first column and the individual detector spectra in subsequent columns. It also adds the background in another column and the standard deviation from counting statistics in another column. If there is more than one detector the sum is also plotted. If there is no valid energy calibration the first column contains the channel number with the first channel being zero. In this case the detectors are not summed and the background may not be very good.

### Argument list:

* Spectrum file (read)
* Plot file (overwritten)

## Primary

Calculate the primary spectrum.

### Argument list:

* Configuration file (read)
* Plot file (overwritten, optional but the only output of the calculation)

## Calculate

Calculates the full spectrum of the first standard in the standards input file.

### Argument list:

* Configuration file (read)
* Standards input file (read)
* Plot file (overwritten, optional but only output of the calculation)

## Compare

Compare a measured spectrum with a full spectrum calculation.

### Argument list:

* Configuration file (read)
* Standards input file (read)
* Spectrum file (read)
* Plot file (required, overwritten)

## Fit

Fit one standard and produce a plot file. Note that the element fit control list is required but does not need any entries. It can be an empty quoted string.

### Argument list:

* Configuration file (read)
* Standards input file (read)
* Element fit control list (required)
* Plot file (required, overwritten)

## Sum

Sum a sequence of spectra and produce a sum spectrum and a maximum value spectrum. See the map action for a description of spectrum file name sequencing. See the bulk sum action for a description of the sum and max value processing. The sum spectrum is quantified after it is computed.

### Argument list:

* Configuration file (read)
* Calibration file (read)
* Spectrum file, first in the sequence (read)
* Element fit control list (required)
* Plot file (required, overwritten)

## Calibrate

Read the standards input file and produce a calibration file of element calibration factors. The calibration file will not be accessed if errors in the input are detected.

### Argument list:

* Configuration file (read)
* Standards input file (read)
* Calibration file (overwritten if no errors)
* Element fit control list (optional but must be present if any options)

## Quantify

Quantify the amount of the given list of elements in an unknown sample from a measured XRF spectrum.

### Argument list:

* Configuration file (read)
* Calibration file (read)
* Spectrum file (read)
* Element fit control list (required)
* Plot file (required, overwritten)

## Map

Reads a sequence of spectra and produces a map file. The map file is a comma-separated values file with the location of each spectrum and the percents of each element in columns, with each spectrum on a separate line. No plot file is produced. The final columns may contain some diagnostic information, such as the goodness-of-fit parameter (reduced chi squared), to help distinguish any poorly-quantified spectra.

### Argument list:

* Configuration file (read)
* Calibration file (read)
* Spectrum file, first in the sequence (read)
* Element fit control list (required)
* Map file (required, overwritten)
