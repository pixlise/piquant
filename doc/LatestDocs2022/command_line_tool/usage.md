---
id: usage
title: Command Line Usage
---

 The actual calculations as well as reading and writing all the files are handled in a command line tool written in C++. This command line tool can be called directly from a terminal window, included in scripts and batch files, and called as a subprocess from other applications. (This is how the GUI accesses it.) Each action listed above has a corresponding sub-command in the argument list.

The text output that is displayed in the text output window of the main window for the GUI is the same text that is produced by the command line tool. That text is written to standard output and will show up in a terminal window. Including a terminal output file will cause this text to go to the file instead. 