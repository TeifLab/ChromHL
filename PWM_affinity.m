% This function computes the affinity across a sequence using the TRAP
% algorithm of Roider et al., 2007
%
% This is a direct conversion of the R script in the tRap package downloaded
% from http://trap.molgen.mpg.de
% with one difference - if the sequence is N, it computes the energy mismatch to be
% the expected value of the energy mismatch given the base GC and AT content
% rather than 0 in the original algorithm
%
% So credit to the original authors - converted to MATLAB by me, GJT, March 2017

function [a_fwd,a_rev] = PWM_affinity(varargin)

% inputs:
% - Seq_DNA DNA sequence [encoded as A=1;C=2;G=3;T=4]; 
% - K0 (defaults to 1e9)
% - Lambda (defaults to 1.5)
% - pwm matrix (defaults to the text file CTCF_matrix_Orlov.txt)
% - gc_constant (defaults to 0.42)
% - pseudo_count (defaults to 1)

% a change of pseudo_count would allow this to be used for PSSM matrices 
% (as a pseudo_count of 1 would be greater than all elements in the matrix)

R0 = 1e9;
Lambda = 1.5;
gc_content = 0.42; % defaults for background ACGT=[0.28,0.21,0.21,0.28]
pseudo_count = 1;

pwm_matrix = 'CTCF_matrix_Orlov.txt';

switch nargin
    case 6
        Seq_DNA = varargin{1};
        R0 = varargin{2};
        Lambda = varargin{3};
        pwm_matrix = varargin{4};
        gc_content = varargin{5};
        pseudo_count = varargin{6};
    case 5
        Seq_DNA = varargin{1};
        R0 = varargin{2};
        Lambda = varargin{3};
        pwm_matrix = varargin{4};
        gc_content = varargin{5};
    case 4
        Seq_DNA = varargin{1};
        R0 = varargin{2};
        Lambda = varargin{3};
        pwm_matrix = varargin{4};
    case 3
        Seq_DNA = varargin{1};
        R0 = varargin{2};
        Lambda = varargin{3};
    case 2
        Seq_DNA = varargin{1};
        R0 = varargin{2};
    case 1
        Seq_DNA = varargin{1};
    case 0
        error('Need to specify at least one argument: the DNA sequence\n');
end

% outputs:
% - a_pwm - forward affinities (per base)
% - a_rev - reverse affinities (per base)
  
Seq_length = length(Seq_DNA);

% import PWM matrix

pwm = importdata(pwm_matrix);
pwm = pwm.data;

motif_length = size(pwm,2);
    
% If we don't have an R0, (denoted by NaN when passed to the function)
% then base it on the equation in the above paper
if isnan(R0)
    R0 = exp(0.584 * motif_length - 5.66);
end
  
% Compute AT content
at_content = 1 - gc_content;
    
% add in the pseudo_count
pwm = pwm + pseudo_count;
    
% This will hold a log transform of the energy-mismatch
% It is all coded in logs to ensure numerical stability
% and reduce inaccuracy
transformed = zeros(4,motif_length);
    
for i=1:motif_length
    max_AT = max(pwm([1 4],i));
    max_GC = max(pwm([2 3],i));
    if (max_AT > max_GC)
        transformed(:,i)=[...
            log(max_AT / pwm(1,i) )/ Lambda, ...
            log((max_AT / at_content) * (gc_content / pwm(2,i)))/ Lambda, ...
            log((max_AT / at_content) * (gc_content / pwm(3,i)))/ Lambda, ...
            log(max_AT / pwm(4,i) )/ Lambda];
    else
        transformed(:,i)=[...
            log((max_GC / gc_content) * (at_content / pwm(1,i)))/ Lambda, ...
            log(max_GC / pwm(2,i) )/ Lambda, ...
            log(max_GC / pwm(3,i) )/ Lambda, ...
            log((max_GC / gc_content) * (at_content / pwm(4,i)))/ Lambda ];
     end
     if (max_AT == max_GC)
         transformed(:,i)= log(max_AT ./ pwm(:,i)) / Lambda;
     end
end
    
% Calculate the reverse-complement log-transformed values
complement = zeros(4,motif_length);
   
for m=1:motif_length
    complement(:,motif_length-m+1) = transformed(4:-1:1,m);
end

% Now calculate the energy mismatches along the sequence
% starting at the left end
n=1;

% These will hold the binding constants (K) and affinities
% (a) of each motif_length-sized window 
K_fwd = zeros (1,Seq_length-motif_length+1);
K_rev = zeros (1,Seq_length-motif_length+1);
a_fwd = zeros (1,Seq_length-motif_length+1);
a_rev = zeros (1,Seq_length-motif_length+1);
    
prob_base = 0.5*[at_content gc_content gc_content at_content];

% Work along the sequence
while (n <= Seq_length-motif_length+1)
    % These variables hold the energy mismatch between sequence and motif
    dE_forward = 0;
    dE_compl = 0;
    
    % work along the motif starting at point n
    for m=1:motif_length
        base = Seq_DNA(n+m-1);
        if (base<5) % if A,C, G or T
            dE_forward = dE_forward + transformed(base,m);
            dE_compl = dE_compl + complement(base,m);
        else % if not A, C, G, T, then use expected value of energy mismatch
            dE_forward = dE_forward + sum(transformed(:,m)'.*prob_base);
            dE_compl = dE_compl + sum(complement(:,m)'.*prob_base);
        end
    end 

    % Calculate binding constant at this point
    K_fwd(n) = R0 * exp (-1*dE_forward);
    K_rev(n) = R0 * exp (-1*dE_compl);
    
    % Calculate affinity by K/(1+K)    
    a_fwd(n) = K_fwd(n) / (1+K_fwd(n));
    a_rev(n) = K_rev(n) / (1+K_rev(n));
        
    n=n+1;
        
end   

end

