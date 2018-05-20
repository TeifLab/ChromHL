% This script runs through sequences in variable Seq_DNA and outputs a MATLAB file
% containing the binding maps, and a text file for each individual profile
%
% Input:
% Sequences (or set of sequences)
% Heterochromatin profile for this set of sequences

noSequences = length(Seq_DNA);
profiles = cell(1,noSequences);

% For each sequence calculate binding map

for i=1:noSequences

    % Set up variables
    seq_to_calculate = Seq_DNA(i).sequence;

    % Mark the centre of the sequence, for correctly centering the lattice by flipping 
    % the value of the sequence in the centre
    seq_to_calculate(Seq_DNA(i).centre_point) = -seq_to_calculate(Seq_DNA(i).centre_point);
    
    % Set the heterochromatin marks for this sequence
    heterochromatin_marks = Seq_DNA(i).het_marks;
    
    % Set text output file name and remove ":" 
    file_out = [Seq_DNA(i).header '.txt'];
    file_out = strrep(file_out, ':', '.');

    % Display output filename
    disp(file_out);
    
    % Calculate the binding map using the variables above
    profiles{i} = CalculateMap(seq_to_calculate,heterochromatin_marks,file_out);

end

% Write profiles to output file with unique identifier computed from
% last sequence in the collection
zz = strrep(Seq_DNA(i).header,':','.');

save(['output-' zz '.mat'],'profiles');

% zz is returned to multiple_run_mcd_20.m, and is used there to rename the output file
