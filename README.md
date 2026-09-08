# Iradina_CEA_Code


## What is it ?
This is the latest source code of Iradina in its CEA's version. This version may be slightly ahead of the code compiled and ditributed in the GUI which in in sourceforge (https://sourceforge.net/projects/iradina/files/Iradina_CEA_GRAPHICAL_USER_INTERFACE/).


##Special features
This version implements 3 different electronic stopping power

stopping_model=O, 1 or 2  with 0 = SRIM, #1 =DPASS ; 2=table.


## Compilation

Details about compilation are in V136/Compilation_Guide

##input file
See the test in TSiC.
Main variables of Config.in are explained in the file.

##electronic stopping
elstop.in is read with stopping_model=2. The first line indicates the number of subsequent lines of data (N). The N lines then give in energy in eV and stopping power in ev/nm.

stopping_model=0 is SRIM 2008.

stopping_model=1 is DPASS formula.

