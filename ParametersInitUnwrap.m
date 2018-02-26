function [] = ParametersInitUnwrap(DNA_seq,heterochromatin_marks)

% input: DNA_seq (per base)
%        heterochromatin_marks (per_base)

global Lpolymer fNumberOfLigands  eNumberOfChromatinStates
global m % ligand lengthes: integer, dimension (1:10)
global mm % (the same as m(g) but real, not integer: double precision, dimension (1:10)
global nMaxGap % max interaction lengthes : integer, dimension (1:10)
global w %cooperativity parameters : double precision w(0:10000,0:10,0:10, 0:5)
global KKK % double precision(10000,10,147,5)
global Unwrap % double precision Unwrap(0:200,1:10)
global c0 % c0(10)
global s sigma % s(1:10000,1:5), sigma(1:5,1:5)
global noGaps

fNumberOfLigands=1; % HP1 ~~and~~CTCF
eNumberOfChromatinStates=3; % heterochromatin, euchromatin and CTCF binding
%Lpolymer=1500;

CTCF_profile = ReadBindingConstants(DNA_seq); % build CTCF affinity curve

het_marks = GenerateHetMarks(DNA_seq); % build het marks (based on PAX3 and PAX9 binding)

%het_marks = [log10(het_marks)>=-7.4]; % number found from plotting TRAP scores (lab book 058)
het_marks = [log10(het_marks)>=-4.5]; % number found from plotting TRAP scores (lab book 058)
%Solid threshold at 10^-1

%het_marks = 2*(log10(het_marks)+2-1) % linear threshold het_marks>10^-1 to het_marks 10^-0.5

het_marks = max(min(het_marks,1),0);



%het_marks = GenerateHetMarksSHIN(DNA_seq); % Use SHIN motif to get initiation sites
%het_marks = (het_marks>(150/160)); % up to 10 mismatches

%het_marks = (het_profile)./(0.5); % now between 0-1
%het_marks = [het_marks>0.9999]; % only keep highest het_marks

%CTCF_scales = (CTCF_profile)/(10^-2); % now scaled so max is 1
%CTCF_scales = [CTCF_scales>0.1];%.*CTCF_scales; % temporary removal of CTCF binding sites

CTCF_scales = [CTCF_profile>1e-5]; % sharp transition to binding sites for affinity 10^-5
%CTCF_scales = CTCF_scales*0; % remove CTCF binding 

Lpolymer = length(CTCF_scales);

s = zeros(Lpolymer,5);
sigma = zeros(5,5);
Unwrap = zeros(201,10); % note that first component was zero-based in F90 code
KKK = zeros(Lpolymer,10,1,5);

s_ratio = 0.2;%9784;
sigma_12 = 1;
cK1 = 0.0132;
w11 = 100;

% parameters for chromatin states: s

if s_ratio>1
    s1=1;
    s2=1/s_ratio;
else
    s1=s_ratio;
    s2=1;
end

s(:,1)=s1;
s(:,2)=s2;

s(:,1)=s_ratio+(1-s_ratio)*het_marks';  % het_marks==1 determines heterochromatin
s(:,2)=s(:,2).*(1-het_marks)';
 
s(:,1)=s(:,1).*(1-CTCF_scales)';
s(:,2)=s(:,2).*(1-CTCF_scales)';
s(:,3)=CTCF_scales;


% sigma
sigma(:,:)=1.;
sigma(1,2)=sigma_12;
sigma(2,1)=sigma_12;

%cost of switching to CTCF binding from either state - is the cost higher or lower than het-eu transition

sigma(1,3)=1;
sigma(3,1)=1;
sigma(2,3)=1;
sigma(3,2)=1;

% sigma(1:2,3)=.000007;
% sigma(3,1:2)=.000007;

% First ligand (HP1)
g=1;
m(g)=1;
mm(g)=real(m(g));
Unwrap(:,g)=0.; %This line should be always equal to zero
Unwrap(1,g)=1.; %This line should be always equal to one

% concentration
c0(g)=1e-7;

% binding constants - preferential binding to state 1 (heterochromatin),
% some binding to state 2 (euchromatin), little or no binding to state 3
% (CTCF)
KKK(:,g,:,1)=cK1/c0(g);
KKK(:,g,:,2)=cK1/c0(g)/100;
KKK(:,g,:,3)=0;


% Cooperativity constants (note indices incremented by 1

w(:,:,:,:)=1.;
w(1,2,2,2:3)=w11;


nMaxGap(1:fNumberOfLigands)=0;
noGaps=false;

end  % subroutine ParametersInitUnwrap
