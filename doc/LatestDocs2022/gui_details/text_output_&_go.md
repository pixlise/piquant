---
id: text_output_&_go
title: Text Output and the GO! Button
---

 Text output is displayed in the text output box at the bottom of the main window. For some actions this is the main information displayed, and for some it is ancillary information and the main output is the plot. The text box is scrollable and editable, although the details of how to edit may vary depending on the operating system.

If a log file is specified, the content of this window is appended to the log file.

None of the inputs are checked or used until the Go! button is pressed. The inputs are then checked to see if all required inputs have been specified. If not, a dialog is displayed showing which input is missing. Once all inputs are verified as present (although not interpreted), the model is invoked to read the input files, perform the requested processing, and write the outputs. Once the model is finished the plot (if any) is displayed in a separate window. The calculations may take several minutes depending on how many channels are in the spectra, how many elements are being quantified, and how many spectra (standards in calibrate or unknowns in map) are being processed. The text output always has the execution time at the end, so future execution times can be estimated. The time increases roughly linearly with the number of channels and spectra and as the square of the number of elements (because of secondary fluorescence cross-calculations). Using the element fit control list can help reduce the execution time as well as improve the fit.

The graphical user interface will be unresponsive after the Go! button is depressed until the action completes and the new text appears in the window. The only way to interrupt the action is to terminate the PIQUANT_CommandLine process using tools provided by the operating system. Note that for large maps the action may take several hours and will give no indication until it completes. Every attempt is made to find errors immediately rather than after a long wait. 