# Iradina_CEA_Code


## What is it ?
This is the latest source code of Iradina in its CEA's version. This version may be slightly ahead of the code compiled and ditributed in the GUI which in in sourceforge (https://sourceforge.net/projects/iradina/files/Iradina_CEA_GRAPHICAL_USER_INTERFACE/).

This code originates in a development by the University of Jena. It has been pursued by French CEA and is GPL licensed https://www.gnu.org/licenses/gpl.html


##Special features
This version implements 3 different electronic stopping power

stopping_model=O, 1 or 2  with 0 = SRIM, #1 =DPASS ; 2=table.


## Compilation

Details about compilation are in V136/Compilation_Guide

##input files and test
###Au_SiC
This test (gold in SiC) has been originally built with GUI. It contains multiple input files. It launches with :

\>iradina -p 9 -data ../data -c ./Configuration.in

Note that the data directory is in AuSiC/../data

###test_UO2
This test uses a single input file (Config.in). It lanches with the command:

\>iradina

In this example the data directory is in the current directory.
Main variables of Config.in are explained in the file.


##electronic stopping
###stopping_model=0 is SRIM 2008.

###stopping_model=1 is DPASS formula (see NIMB 460, 19 (2019). The file DPASS_DB.dat should be in the data directory (see below)

###elstop.in is read with stopping_model=2.
The file is built as sequence of blocks, one block per element involved in the calculation, with the target elements then the projectile element, i.e. for Au in SiC, the first blcok is for Si, then C then Au.

For each block, the first line indicates the number of subsequent lines of data (N). The subsequent N lines then give in energy in eV and stopping power in ev/nm for this element in the present target.

This file has to be built or obtained from outside iradina. It must be in the working directory, not in data/ (as it is very case specific).

