     TALYS-1.6 (Version: December 23, 2013)

 Copyright (C) 2013  A.J. Koning, S. Hilaire and S. Goriely        
                         NRG          CEA              ULB

The TALYS package
-----------------

In what follows we assume TALYS will be installed on a Unix/Linux 
operating system. The total TALYS package is in the talys/ directory and 
contains the following directories and files:

- README: this file

- talys.setup is a script that takes care of the installation.

- source/ contains the source code of TALYS: 274 Fortran subroutines, and the 
  file talys.cmb, which contains all variable declarations and common blocks.
  This includes the file ecis06t.f. This is basically Jacques Raynal's code 
  ECIS-06, which we have transformed into a subroutine and slightly modified 
  to generate extra output that is not given by the original version of ECIS.

- structure/ contains the nuclear structure database in various subdirectories. 

- doc/ contains the documentation: this manual in postscript and pdf format and 
  the description of ECIS-06.

- samples/ contains the input and output files of the sample cases.

In total, you will need about 4 Gb of free disk space to install TALYS.

Installation
------------

The installation of TALYS is straightforward.
For a Unix/Linux system, the installation is expected to be handled by the 
talys.setup script. This can be activated by

- edit talys.setup and set the first two variables: the name of your compiler 
  and the place where you want to store the TALYS executable.

- talys.setup

If this does not work for some reason, we here provide the necessary steps to 
do the installation manually. For a Unix/Linux system, the following steps 
should be taken:

- chmod -R u+rwX talys 
  (to ensure that all directories and files have read and write permission and 
  the directories have execute permission. 

- cd talys/source

- Ensure that TALYS can read the nuclear structure database. This is done
  in subroutine machine.f. If talys.setup has not already replaced the path 
  name in machine.f, do it yourself. We think this is the only Unix/Linux 
  machine dependence of TALYS. Apart from a few trivial warning messages for 
  ecis06t.f, we expect no complaints from the compiler. 

- f95 -c *.f

- f95 *.o -o talys

- mv talys /bin 
  (assuming you have a /bin directory which contains all executables that can 
  be called from any working directory)

After you or talys.setup has completed this, type

- rehash
  (to update your table with commands)

Verification
------------

- cd samples

- verify

Under Linux/Unix, this should run all sample cases (about 1 hour on
a fast PC).

Should you encounter error messages, upon running TALYS, like 'killed' or
'segmentation fault', then probably the memory of your processor is not large
enough (i.e. smaller than 256 Mb). Edit talys.cmb and change the value of 
memorypar.

Your own calculations
---------------------

- talys < input > output

where you can make your own input file starting from the many sample cases
we provide.
