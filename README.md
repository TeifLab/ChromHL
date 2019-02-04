# ChromHL

Code for computation of chromatin microdomains based on sequence.

## Overview:

Code written in combination of R and MATLAB for calculation of chromatin microdomains based upon either transcription factor binding or repeat recognition.

#### Main loop runs as follows:
 
1. Initialise heterochromatin locations (either by calculating binding affinity or by location of relevant genomic feature, such as repeat)
2. Set up CTCF binding sites by calculating binding affinity
3. Compute states based on initiation sites and CTCF binding sites, based on biophysical parameters such as the far-field statistical weights for each chromatin state, the statistical weight of forming a boundary between two chromatin states, the cooperativity of the associated protein HP1 and the binding constant of HP1 to each chromatin state

#### Required inputs are:

1. Sequence of region to simulate.
2. Location of relevant repeats (as part of a BED-like file), if repeats are being used as initiation sites.

#### Outputs are:

1. Concentration (in units of initial HP1 concentration) of HP1 binding
2. Probability of each chromatin state per lattice unit.

#### Note that:

1. Lattice units are often taken to be approximately equal to the NRL for the region of study, so ~179-189bp.
2. Algorithm for TF-binding recognition uses a sliding geometric-average window of 501bp, as windows smaller than this are too noisy in terms of predicting heterochromatin from TF-binding

Other things to note are that code is based on an original FORTRAN version from @epigenereg, converted into MATLAB by @geejaytee, with some scripts for analysis in R by @geejaytee.
