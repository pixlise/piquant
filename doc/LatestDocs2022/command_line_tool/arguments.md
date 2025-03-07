---
id: arguments
title: Command Line Arguments
---

 The sub-command determines the rest of the argument list and how it is interpreted. The preceding sections give the expected argument list for each sub-command. An error results if the argument list does not match what is expected for the sub-command.

All file arguments are full path. The system file open is used so its rules apply for file paths. See the section on file formats for a description of each type of file.

The full path argument for the plot file is required when it is used because there is no provision for a temporary file and the argument list is order-dependent. It can be omitted from some sub-commands, where it is optional, if there are no arguments following it. Note that there is no display of the plot, only writing the CSV plot file.

The element fit control list must be enclosed in quotes if it contains any blanks or horizontal tabs. Once the argument is passed to the command line tool it is parsed just as if it came from the graphical user interface, so its use is as described above. Note that it is optional for some actions in the GUI but must always be present in the command line arguments if it is listed under a sub-command (although it can sometimes be an empty quoted string).

A full path argument for the terminal output file can always be appended to the argument list. It must follow all other, required arguments and must precede any options. The terminal output file can be used as a log file because it includes a header and time stamp. It always appends and most operating systems will create a file opened for append if one does not already exist. 