% Driver script for running map

% Input:
% Sequences (or set of sequences)
% Heterochromatin profile for this set of sequences

noSequences = length(Seq_DNA);
profiles = cell(1,noSequences);

% For each sequence
% Calculate map

for i=1:noSequences
    seq_to_calculate = Seq_DNA(i).sequence;
    seq_to_calculate(Seq_DNA(i).centre_point) = -seq_to_calculate(Seq_DNA(i).centre_point); % flip centre point of sequence
    heterochromatin_marks = Seq_DNA(i).het_marks;
    
    file_out = [Seq_DNA(i).header '.txt'];

    file_out = strrep(file_out, ':', '.');

    disp(file_out);
    
    profiles{i} = CalculateMap(seq_to_calculate,heterochromatin_marks,file_out);
end

% write output to file

zz = strrep(Seq_DNA(i).header,':','.');

save(['output-' zz '.mat'],'profiles');

