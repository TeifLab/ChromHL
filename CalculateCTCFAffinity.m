% This function reads in the DNA_seq and calculates the binding affinity
% map for CTCF using the PWM method, divided into lattice_size bp 
% windows

function [CTCF_profile] = ReadBindingConstants(DNA_seq)

% Set up lattice size and two constants based on lattice
lattice_size=161;
half_size_low=floor(lattice_size/2);
half_size_high=half_size_low+1;

% Find where the DNA_seq variable is negative, as this is where the centre
% of the sequence has been marked
numbers = (1:length(DNA_seq));
centre_point = numbers(DNA_seq<0); % centre point is marked by negative seq value

% Now we don't need the negative element, so flip sign back to positive
DNA_seq = abs(DNA_seq);

% Calculate binding affinities in forward and reverse directions
[K_pwm,K_rev] = PWM_affinity(DNA_seq,NaN,0.7,'CTCF_matrix_Orlov.txt');

% Smooth out the affinities across the length of the binding site
% - this is an approximation to the full binding calculation of CTCF

% first find the length of the motif
pwm = importdata('CTCF_matrix_Orlov.txt');
pwm = pwm.data;
motif_length = size(pwm,2);
    
K_fwd=[K_pwm zeros(1,motif_length-1)];
K_rev=[K_rev zeros(1,motif_length-1)];

% build shifts of the affinities across the motifs        
matrix_fwd = K_fwd;
matrix_rev = K_rev;
        
for j=1:motif_length-1
    matrix_fwd=[matrix_fwd; circshift(K_fwd,j)];
    matrix_rev=[matrix_rev; circshift(K_rev,j)];
end

% take the maximum "vertically" to build the approximation to the full
% binding map calculation
K_fwd=max(matrix_fwd,[],1);
K_rev=max(matrix_rev,[],1);

% CTCF can bind to either strand
K_sum = K_fwd+K_rev.*(1-K_fwd);

% Now convert the CTCF binding at bp level to binding at lattice unit level

% First count the number of units either side of the centre point
left_nucleosomes = ceil((centre_point-half_size_high)/lattice_size);
right_nucleosomes = ceil((length(K_sum)-centre_point-half_size_low)/lattice_size);
% then the number of nucleosomes in string: left_nucleosomes + 1 + right_nucleosomes
% This ensures that the sequence is completely covered by lattice units

% Add zeros to either side of the affinity curve
K_sum = [zeros(1,left_nucleosomes*lattice_size-(centre_point-half_size_high)) ...
    K_sum zeros(1,right_nucleosomes*lattice_size-(length(K_sum)-centre_point-half_size_low))];

% Reshape so that each column of K_sum is the affinities in each lattice unit
K_sum = reshape(K_sum,[lattice_size left_nucleosomes+1+right_nucleosomes]);

% Sum along the columns to give the affinities per lattice unit
CTCF_profile = sum(K_sum,1);

end
