---
id: file_selection
title: File Selection Inputs
---

 Most of the inputs are in the form of files prepared in advance and then selected using the main window of the graphical user interface. File formats are described in Section 7. Many of the files are produced automatically, such as the spectrum files that are produced during data collection. Push the select button to bring up a system file dialog to point to the location of the file. The full path to the file will be placed into the text line under the select button. This text can be edited to change the path if desired or a path can be entered manually (or via copy and paste). The dialogs will show only the correct file type for the input being selected.

The text line for the element fit control list is different in that it can only be entered manually (or via copy and paste). It is direct text sent to the PIQUANT calculation engine and is not a file path. The text consists of a list of element symbols or atomic numbers and some single-letter instructions. The actions that use this input vary in how it is interpreted.

The plot file selection is optional. If a plot will be produced and a plot file is not specified, the plot information will be written to a temporary file and the plot produced from this file. The temporary file will be deleted when the plot is completed.

The log file is always available. If a log file name is specified, the text output that appears in the output text window above the Go! button will be appended to the log file after a date and time stamp. A new log file can be created or an existing file chosen, depending on the choice beside the select button. In either case the log entries will always be appended, never overwritten. Note that the log file is written after the text is displayed in the window, so nothing new will be written to the log file until the action completes. This may take a while for big maps. 