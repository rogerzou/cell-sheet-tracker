cell-sheet-tracker
================

### Citing our work
This source code accompanies the paper below. Please cite if you use our code for your research. All code is written in the MATLAB and Java programming languages.

R. S. Zou and C. Tomasi. Deformable graph model for tracking epithelial cell sheets in fluorescence microscopy. *IEEE Transactions on Medical Imaging*, 35(7):1625-1635, 2016.

### Software requirements
- MATLAB R2015b (wth Image Processing and Statistics Toolboxes)
- Java Development Kit (JDK). Use an older JDK version than the Java version that comes with your MATLAB distribution (found using `version -java`)
- SIFT flow <http://people.csail.mit.edu/celiu/SIFTflow/> (optional, see Section II-F of the paper)

### Installation and setup
1. Compile the source files in `src/java/` into a new folder `bin/` with the appropriate JDK (see `compile.sh`)
2. If using SIFT flow, download from the link above and build according to their instructions. Place the SIFT flow folder into a new folder `lib/`.
3. In MATLAB, run `setup`.

### OS X Specific Notes
- Running `/usr/libexec/java_home -V` on command line lists versions of java installed on your system. To use to a specific version, run `export JAVA_HOME="path/to/jdk"`
- To compile java files, try running `bash compile.sh` on command line in the home directory.

### Demonstrations
The script `RUNME` demonstrates our algorithm for a subset of the DDC1 data set. If you decide to not use SIFT flow, change the `siftflow` flag in line 26 of `RUNME` to `false`.
