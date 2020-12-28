# Program flow

## Introduction

This program flow is relative to running in standalone mode, running the script ````multiple_run_all.m```` inside an interative MATLAB session. For running on a cluster, submitting jobs via ````job_submit.sh````, then the flow is the same, starting from when ````multiple_run_chromhl.m```` is run.

## Initialisation

  - The script ````multiple_run_<type>.m```` reads in all FASTA sequences using ````Read_FASTA_all.m```` and runs ````Driver.m```` to perform the main loop
  
## Main loop
 
  - For each sequence provided in, ````Driver.m```` computes the binding maps using ````CalculateMap.m```` and outputs a MATLAB file containing the maps, and a text file for each individual profile
  - ````CalculateMap.m```` runs in four parts:
      1. ````SetDefaults.m```` sets parameter defaults
      2. ````ParametersInitMicrodomain.m```` initialises the statistical weights, the initiation sites and the CTCF-occupied lattice units
        
          - ````CalculateCTCFAffinity.m```` applies the TRAP algorithm (via ````PWM_affinity.m````) on the sequence level, then takes the maximum affinity per lattice unit to give a lattice-based affinity.
        
          - ````CalculatePAX39Affinity.m```` applies the TRAP algorithm (via ````PWM_affinity.m````) for PAX3 and PAX9 transcription factors, then smooths the affinity using a geometric mean over a 501-bp centred windw
        Once affinities are computed per lattice unit, the lattice unit is set as heterochomatin initiation site if the lattice PAX3/9 affinity is above a threshold or as CTCF-bound if the lattice CTCF affinity is above a threshold.
      3. ````ConstantsInitMicrodomain.m```` fixes the transfer matrix dimension and enumerates lattice unit states
      4. ````MapOfBindingCalc.m```` computes the binding map using ````PointOfMapOfBindingCalc.m````, running along the lattic
      
          - ````PointOfMapOfBindingCalc.m```` computes the HP1 binding probability, and the probability of each chromatin state at that point in the lattice
              - Initialising the transfer matrix and its derivatives using ````MatrixInitMicrodomain.m````, and computing the relevant state variables at that point in the lattice, based on the transfer matrix at all other points in the lattice

## Outputs

  - The outputs are a MATLAB MAT file containing the binding maps for each sequence in the input MAT file, and a fixed-format text file with the binding map per lattice unit
  
## Expected running time

  - For the four 8000-bp sequences provided as demonstration, the calculation of all maps takes approximately 5 seconds on an Intel Core i7-8550U CPU with 8 GB RAM.
