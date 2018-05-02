function [CTCF_profile] = ReadBindingConstants(DNA_seq)

% function uses read in DNA_seq to calculate binding constants for CTCF
% using PWM method
% Note also that it is divided into lattice_size-bp windows
% so CTCF_profile = [prod(1,187) K] ^ (1/lattice_size) as geometric mean

%%% need centre point - centre_point

% call PWM method for base-pair level binding constants

lattice_size=161;
half_size_low=floor(lattice_size/2);
half_size_high=half_size_low+1;

numbers = (1:length(DNA_seq));

centre_point = numbers(DNA_seq<0); % centre point is marked by negative seq value

DNA_seq = abs(DNA_seq);

[K_pwm,K_rev] = PWM_affinity(DNA_seq,NaN,0.7,'CTCF_matrix_Orlov.txt');

% smooth out the affinities across binding sites

pwm = importdata('CTCF_matrix_Orlov.txt');
	
pwm = pwm.data;

motif_length = size(pwm,2);
    

K_fwd=[K_pwm zeros(1,motif_length-1)];
K_rev=[K_rev zeros(1,motif_length-1)];
        
%        K_fwd(isinf(K_fwd))=NaN;
%        K_rev(isinf(K_rev))=NaN;
        
        matrix_fwd = K_fwd;
        matrix_rev = K_rev;
        
        for j=1:motif_length-1
            matrix_fwd=[matrix_fwd; circshift(K_fwd,j)];
            matrix_rev=[matrix_rev; circshift(K_rev,j)];
        end
        
        K_fwd=max(matrix_fwd,[],1);
        K_rev=max(matrix_rev,[],1);

K_sum = K_fwd+K_rev; % can bind to either strand

%K_sum = K_sum+flip(K_sum);

%width=20*147;

%K_sum = 1-exp(-(-20000:20000).^2/(2*width^2));


left_nucleosomes = ceil((centre_point-half_size_high)/lattice_size);
right_nucleosomes = ceil((length(K_sum)-centre_point-half_size_low)/lattice_size);

% number of nucleosomes in string: left_nucleosomes + 1 + right_nucleosomes
% 
% K_sum = [nan(1,left_nucleosomes*lattice_size-(centre_point-half_size_high)) ...
%   K_sum nan(1,right_nucleosomes*lattice_size-(length(K_sum)-centre_point-half_size_low))];
% 
% 
% K_sum = reshape(K_sum,[lattice_size left_nucleosomes+1+right_nucleosomes]);
% 
% CTCF_profile = exp(nanmean(log(K_sum),1));

% number of nucleosomes in string: left_nucleosomes + 1 + right_nucleosomes

K_sum = [zeros(1,left_nucleosomes*lattice_size-(centre_point-half_size_high)) ...
    K_sum zeros(1,right_nucleosomes*lattice_size-(length(K_sum)-centre_point-half_size_low))];


K_sum = reshape(K_sum,[lattice_size left_nucleosomes+1+right_nucleosomes]);

CTCF_profile = sum(K_sum,1); % As in TRAP calculation affinity per lattice_size-bp region is sum of affinities for each bp in region


end
