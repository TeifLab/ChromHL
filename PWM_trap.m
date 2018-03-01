    %%%
    % all definitions as in Roider et al., 2007 (Vingron)
    % max_i(i) - maximal element in the column
    % dE(i, base) - energy change for base pair with respect to the strongest site
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    function [a_fwd,a_rev] = PWM_affinity(varargin)

    % inputs:
    % - Seq_DNA DNA sequence [encoded as A=1;C=2;G=3;T=4];
    % - K0 (defaults to 1e-9)
    % - Lambda (defaults to 0.7)

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
    
    if isnan(R0)
        R0 = exp(0.584 * motif_length - 5.66);
    end
  
    at_content = 1 - gc_content;
    
    pwm = pwm + pseudo_count;
    
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
    
    complement = zeros(4,motif_length);
    
    for m=1:motif_length
        complement(:,motif_length-m+1) = transformed(4:-1:1,m);
    end
    
    %P_combined = 0;
    %P_uncorrected = 0;
    
    n=1;
    
    K_fwd = zeros (1,Seq_length-motif_length+1);
    K_rev = zeros (1,Seq_length-motif_length+1);
    a_fwd = zeros (1,Seq_length-motif_length+1);
    a_rev = zeros (1,Seq_length-motif_length+1);
    
    prob_base = 0.5*[at_content gc_content gc_content at_content];
      
    while (n <= Seq_length-motif_length+1)
        dE_forward = 0;
        dE_compl = 0;
        
        for m=1:motif_length
            base = Seq_DNA(n+m-1);
            if (base<5)
                dE_forward = dE_forward + transformed(base,m);
                dE_compl = dE_compl + complement(base,m);
            else
                dE_forward = dE_forward + sum(transformed(:,m)'.*prob_base);
                dE_compl = dE_compl + sum(complement(:,m)'.*prob_base);
            end
        end
        
        K_fwd(n) = R0 * exp (-1*dE_forward);
        K_rev(n) = R0 * exp (-1*dE_compl);
        
        a_fwd(n) = K_fwd(n) / (1+K_fwd(n));
        a_rev(n) = K_rev(n) / (1+K_rev(n));
        
        n=n+1;
        
    end   

%    fid15 = fopen('K.txt','w');

%    for i=1:(Seq_length-motif_length+1)
%        fprintf(fid15,'K_fwd(%d)=%f\n',i,K_fwd(i));
%    end
%    for i=1:(Seq_length-motif_length+1)
%        fprintf(fid15,'K_rev(%d)=%f\n',i,K_rev(i));
%    end
%
%    fclose(fid15);

    end

