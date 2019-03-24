% This function sets up the initial parameters for the simulation
% such as the binding sites for the transcription factors PAX3/9
% the binding sites for CTCF
% and
% sigma(i,j)        : overall energy costs of forming boundaries between states i and j
% s(n,e)            : statistical weights of each state at each lattice unit
% KKK(n,g,l,e)      : binding constants for each protein g to each state (element 4) at each lattice unit n, 
%                   : proteins have binding length 1:max(l)
% w(n+1,g+1,h+1,e+1): cooperativity of binding between proteins g and h at a distance n for state e
%                   : Note that due to MATLAB not using 0-based arrays, each index in w is incremented
%                   : by 1
% Unwrap(i+1,g)     : Chance of unwrapping of protein g, again originally 0-based
%
% Input is the DNA_sequence, and the (given) heterochromatin_marks

function [] = ParametersInitMicrodomain(DNA_seq,heterochromatin_marks)

global Lpolymer fNumberOfLigands  eNumberOfChromatinStates
global m % ligand lengths
global nMaxGap % max interaction lengths 
global w % cooperativity parameters  
global KKK % binding constants
global Unwrap % unwrapping constants
global c0 % concentrations
global s sigma % s for states per lattice unit, sigma for forming boundaries between lattice units
global noGaps

% Change these lines for more binding ligands and more chromatin states
% As it stands, only 1 ligand (HP1) 
% and 3 states (1: heterochromatin, 2: euchromatin, 3: CTCF binding)
fNumberOfLigands=1; 
eNumberOfChromatinStates=3;

% Build CTCF affinity profile
CTCF_profile = CalculateCTCFAffinity(DNA_seq); 

% Build PAX39 affinity profile
PAX39_profile = CalculatePAX39Affinity(DNA_seq); 

% Threshold PAX39 affinity to give positions where heterochromatin initiates
het_marks = [log10(PAX39_profile)>=-4.5];

% If the above changes, such as using a linear het_mark probability, ensure that it lies between 0 and 1
het_marks = max(min(het_marks,1),0);

% Threshold CTCF affinity profile to give CTCF binding
CTCF_scales = [CTCF_profile>1e-5]; % sharp transition to binding sites for affinity 10^-5

% Length of calculation domain
Lpolymer = length(CTCF_scales);

% Initialise all variables
s = zeros(Lpolymer,5);
sigma = zeros(5,5);
Unwrap = zeros(201,10); 
KKK = zeros(Lpolymer,10,1,5);

% Initial values for parameters: 
% s_ratio : statistical weight of forming heterochromatin relative to euchromatin
% sigma_12: exp(-delta * energy cost of forming boundary/ RT), so always <=1
% cK1     : concentration * binding constant of HP1
% w11     : cooperativity of HP1 on neighbouring lattice units
s_ratio = 0.2;
sigma_12 = 1;
cK1 = 0.0132;
w11 = 100;


% Set individual state statistical weights
% for numerical stability, max(s)=1
if s_ratio>1
    s1=1;
    s2=1/s_ratio;
else
    s1=s_ratio;
    s2=1;
end

% Set up "background" statistical weights of states
s(:,1)=s1;
s(:,2)=s2;

% Alter so that s(i,1)=1; s(i,2)=0 if site is a heterochromatin nucleation site
s(:,1)=s_ratio+(1-s_ratio)*het_marks'; % het_marks==1 is a definite nucleation site
s(:,2)=s(:,2).*(1-het_marks)';

% Set up CTCF binding sites. Note that CTCF binding "beats" heterochromatin formation
s(:,1)=s(:,1).*(1-CTCF_scales)';
s(:,2)=s(:,2).*(1-CTCF_scales)';
s(:,3)=CTCF_scales;

% Set up energy barriers between states 1 (het) and 2 (eu) from given sigma_12 parameter
sigma(:,:)=1.;
sigma(1,2)=sigma_12;
sigma(2,1)=sigma_12;

% If there is extra cost to forming a CTCF binding site from either state, alter these lines
sigma(1,3)=1;
sigma(3,1)=1;
sigma(2,3)=1;
sigma(3,2)=1;

% Now set up ligands

% First ligand (HP1)
g=1;
m(g)=1; % length of protein in lattice units

% as there is no unwrapping of HP1 binding
Unwrap(:,g)=0.; 
Unwrap(1,g)=1.; 

% Cncentration of HP1
c0(g)=1e-7;

% Set up binding constants - preferential binding to state 1 (heterochromatin),
% some binding to state 2 (euchromatin), no binding to state 3
% (CTCF)
KKK(:,g,:,1)=cK1/c0(g);
KKK(:,g,:,2)=cK1/c0(g)/100;
KKK(:,g,:,3)=0;


% Cooperativity constants (note indices incremented by 1 for all indices)
w(:,:,:,:)=1.;
w(1,2,2,2:3)=w11; 
% So interpreting the above:
% between protein 1 (HP1) and protein 1 (HP1)
% bound in either state 1 or 2, the cooperativity statistical weight is w11

% As there is no minimum gap between ligands
nMaxGap(1:fNumberOfLigands)=0;
noGaps=false;

end 
