cell-sheet-tracker
================

### Citing our work

This source code accompanies the paper below. Please cite if you use our code for your research. All code is written in the MATLAB and Java programming languages.

Zou R.S. and Tomasi, C. (under review). Deformable Graph Model for Tracking Epithelial Cell Sheets in Fluorescence Microscopy. *Medical Imaging, IEEE Transactions on*

### Software requirements
- MATLAB R2015b (wth Image Processing and Statistics Toolboxes)
- Java Development Kit (JDK). Use an older JDK version than the Java version that comes with your MATLAB distribution (found using `version -java`)
- SIFT flow <http://people.csail.mit.edu/celiu/SIFTflow/> (optional, see Section II-F of the paper)

### Installation and setup
1. Compile the source files in `src/java/` into a new folder `bin/` with the appropriate JDK.
2. If using SIFT flow, download from the link above and build according to their instructions. Place the SIFT flow folder into a new folder `lib/`.
3. In MATLAB, run `setup`.

### Demonstrations
The script `RUNME` demonstrates our algorithm for a subset of the DDC1 data set. If you decide to not use SIFT flow, change the `siftflow` flag in line 26 of `README` to `false`.
