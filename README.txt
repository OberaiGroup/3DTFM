
The code is organized as follows:
- folder Elibrary contains the elements files,
- folder Exe      contains the Makefile and the memory allocation file for the modules,
- folder Source   contains the source code for the optimization,
- folder examples contains a series of examples in 2D and 3D for testing the code and
                           learn how to write input files.

Once the executable is created (cf Exe/README.txt), the code is run by typing:
${PathToTheExecutable}nlace.exe ${PathToInputFile}NameOfInputFile.in

The executable outputs the following text files:
*.vtk: result files for paraview
*.res: restart files for NLACE
*.cvg: check the convergence of the optimization
*.ite: if running L-BFGS-B optimization, this file details the BFGS progess


