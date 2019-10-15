
NLACE has been successfully compiled on Linux and MacosX platforms.

NLACE needs to be linked with several exterior libraries which are detailed in the Makefile.
Beyond the blas and lapack libraries, the code needs an implementation of the pardiso solver.
Pardiso is implemented in the Intel MKL libraries or a stand alone library can be downloaded
from the Basel University's website. Both link and execute properly with NLACE.

As of 2010/04/01, NLACE implements OpenMP parallelization by default. This has an impact on 
the compilers that can be used (the code requires OpenMP3.0 and above):
- for GNU, please use gcc/gfortran4.4.0 and above
The following compilers need to be tested but should work:
- for Intel, please use ifort11.0 and above
- for PGI, please use pgfortran10.4 and above

The Makefile has enough comments to be modified so as to run on various platforms. When
the path are properly written in the Makefile, the executable is created by typing:
make



