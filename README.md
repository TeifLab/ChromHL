# ChromHL [![CC BY-NC 4.0][cc-by-nc-shield]][cc-by-nc]

The mammalian epigenome contains thousands of heterochromatin nanodomains (HNDs) marked by di- and trimethylation of histone H3 at lysine 9, which have a typical size of 3-10 nucleosomes. However, the (epi)genetic determinants of their location and boundaries are only partly understood. Here, we compare four HND types in mouse embryonic stem cells, that are defined by histone methylases SUV39H1/2 or GLP, transcription factor ADNP or chromatin remodeller ATRX. Based on a novel chromatin hierarchical lattice framework termed **ChromHL**, we are able to predict HND maps with singe-nucleotide resolution. We find that HND nucleation can be rationalized by DNA sequence specific protein binding to PAX3/9, ADNP and LINE1 repeats. Depending on type of microdomains, boundaries are determined either by CTCF binding sites or by nucleosome-nucleosome and nucleosome-HP1 interactions. Our new framework allows predicting how patterns of H3K9me2/3 and other chromatin nanodomains are established and changed in processes such as cell differentiation.

![alt text](https://github.com/TeifLab/ChromHL/blob/master/GitHub_Header_v2.png?raw=true)

<p align="center">
  <img src="https://github.com/TeifLab/ChromHL/blob/master/GitHub_Header_v2.png">
</p>

## Introduction

This code is written in combination of R and MATLAB. It allows calculations of the epichromatin landscape for a individual genomic region, as well as bulk calculations for a set of genomic regions defined in the FASTA file.

## System requirements
### For running the model and analysing the outputs
- MATLAB - https://uk.mathworks.com/ (tested on >R2018a); this should also work on Octave (https://www.gnu.org/software/octave) instead of MATLAB, but it was not tested.
- R - https://cran.r-project.org/ (tested on >v2.12)
### For generating sequences to compute over
- bedtools - https://bedtools.readthedocs.io/ (tested on >v2.24)

## Installation instructions
- Clone the github repository.

### For standalone running on a few regions (including the demo)
- From a MATLAB command line, run ````multiple_run_all.m````.
    - This will generate one text file per input FASTA file (four, in the case of the demo) containing the predicted HP1 binding, and the probabilities of being in the three chromatin states (heterochromatin, euchromatin, bound to CTCF).
    
### For a task array job on a SGE
- Make a directory ````individual_regions```` containing the individual FASTA files as ````gene_region_XXXX.fa````.
- Alter the header ````-t```` in ````job_submit.sh```` for the number of tasks (ceiling(number of regions/1000)).
- Change the other submission parameters in the header of the ````job_submit.sh```` script.
- Submit the job script ````job_submit.sh````.
    - The text output files will be generated in a directory ````Text_output```` ready for processing.
    
### For a basic visualisation
- From an R command line, run ````plot_chr_states.R```` inside the directory containing the text outputs.

## Code overview
### Main calculation loop
 
1. Initialise heterochromatin locations (either by calculating binding affinity or by location of relevant genomic feature, such as repeat)

    a. Binding affinity of initiation factors is calculated using a sliding window along the sequence, which is then averaged using a geometric mean across a 501-bp centred window. If the mean affinity in a lattice unit is greather than some fixed threshold, then the lattice unit is initialised as being in heterochromatin state, otherwise it remains as euchromatin.
    
    b. For repeats or other genomic features, the lattice unit where the centre of the feature is located is initialised as being in heterochromatin state
2. Set up CTCF binding sites by calculating binding affinity

    a. As for the initiation factors, the binding affinity for CTCF is calculated using a sliding window, without any further averaging; if the affinity within a lattice unit crosses a fixed threshold, the unit is initialised as the CTCF-bound state.
3. Compute states based on initiation sites and CTCF binding sites, based on biophysical parameters such as the far-field statistical weights for each chromatin state, the statistical weight of forming a boundary between two chromatin states, the cooperativity of the associated protein HP1 and the binding constant of HP1 to each chromatin state

More detailed information on the calculation loop is here: https://github.com/TeifLab/ChromHL/blob/master/Program_Flow.md

### Required inputs

1. The sequences of the region(s) to simulate.
2. The weight matrices for all proteins binding to the sequence that are taken into account.
3. The location of sequence repeats (as part of a BED-like file), if repeats are being used as initiation sites.

### Outputs

1. Concentration (in units of initial HP1 concentration) of HP1 binding.
2. Probability of each chromatin state per lattice unit.

### Technical notes

1. Lattice units are often taken to be approximately equal to the NRL for the region of study, so ~179-189bp. This is set in ````CalculatePAX39Affinity.m```` and ````CalculateCTCFAffinity.m````.
2. The initiation of heterochromatin based on TF-binding is currently in a sliding geometric-average window of 501bp; this can be changed by the user in ````CalculatePAX39Affinity.m````.
3. The code is based on the Fortran version from @epigenereg, converted into MATLAB by @geejaytee, with some scripts for analysis in R by @geejaytee.

## How to cite

Thorn G.J., Clarkson C.T., Rademacher A., Mamayusupova H., Schotta G., Rippe K., Teif V.B. (2020) DNA sequence-dependent formation of heterochromatin nanodomains. bioRxiv 2020.12.20.423673 https://www.biorxiv.org/content/10.1101/2020.12.20.423673v1


## License
This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[![CC BY-NC 4.0][cc-by-nc-image]][cc-by-nc]

[cc-by-nc]: http://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://licensebuttons.net/l/by-nc/4.0/88x31.png
[cc-by-nc-shield]: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg
