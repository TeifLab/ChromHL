% Start the timer

tic;
     
% Read the sequences matching the pattern below (in the [])

% Change this if FASTA sequences are elsewhere

Seq_DNA=Read_FASTA_all(['*.fa']);

% If we have sequences (in the Seq_DNA structure)

if ~isempty(Seq_DNA)

  % Pass to individual binding map calculator

  Driver; % This outputs an output file output-{zz}.mat where {zz} is a variable passed back

end

toc; % Stop timer and output time elapsed - useful for timing runs and optimisation of algorithm








