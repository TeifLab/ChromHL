# ChromHL

The mammalian epigenome contains thousands of heterochromatin nanodomains (HNDs) marked by di- and trimethylation of histone H3 at lysine 9, which have a typical size of 3-10 nucleosomes. However, the (epi)genetic determinants of their location and boundaries are only partly understood. Here, we compare four HND types in mouse embryonic stem cells, that are defined by histone methylases SUV39H1/2 or GLP, transcription factor ADNP or chromatin remodeller ATRX. Based on a novel chromatin hierarchical lattice framework termed **ChromHL**, we are able to predict HND maps with singe-nucleotide resolution. We find that HND nucleation can be rationalized by DNA sequence specific protein binding to PAX3/9, ADNP and LINE1 repeats. Depending on type of microdomains, boundaries are determined either by CTCF binding sites or by nucleosome-nucleosome and nucleosome-HP1 interactions. Our new framework allows predicting how patterns of H3K9me2/3 and other chromatin nanodomains are established and changed in processes such as cell differentiation.

## Overview:

This code is written in combination of R and MATLAB. It allows calculations of the epichromatin landscape for a individual genomic region, as well as bulk calculations for a set of genomic regions edfined in the FASTA file.

#### Main calculation loop:
 
1. Initialise heterochromatin locations (either by calculating binding affinity or by location of relevant genomic feature, such as repeat)
2. Set up CTCF binding sites by calculating binding affinity
3. Compute states based on initiation sites and CTCF binding sites, based on biophysical parameters such as the far-field statistical weights for each chromatin state, the statistical weight of forming a boundary between two chromatin states, the cooperativity of the associated protein HP1 and the binding constant of HP1 to each chromatin state

#### Required inputs:

1. The sequences of the region(s) to simulate.
2. Weight matrices for TFs that are taken into account
3. Location of sequence repeats (as part of a BED-like file), if repeats are being used as initiation sites.

#### Outputs:

1. Concentration (in units of initial HP1 concentration) of HP1 binding
2. Probability of each chromatin state per lattice unit.

#### Technical notes:

1. Lattice units are often taken to be approximately equal to the NRL for the region of study, so ~179-189bp.
2. The initiation of heterochromatin based on TF-binding is currently in a sliding geometric-average window of 501bp; this parameters can be changed by the user
3. Historically, the code is based on the Fortran version from @epigenereg, converted into MATLAB by @geejaytee, with some scripts for analysis in R by @geejaytee.

#### How to cite:

Thorn G.J., Clarkson C.T., Rademacher A., Mamayusupova H., Schotta G., Rippe K., Teif V.B. (2020) DNA sequence-dependent formation of heterochromatin nanodomains. bioRxiv 2020.12.20.423673 https://www.biorxiv.org/content/10.1101/2020.12.20.423673v1

