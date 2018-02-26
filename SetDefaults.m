%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% w(i,j) is the contact interaction between the ligands of type i and j
% Polymer unit types: A,T,G,C,X,Y
% ligand types (or binding mode types): g = 1,..,f
% m(g) is the length of the ligand of type g (in DNA units)
% K(n,g) is the binding constant for type g protein at DNA site [n , n+m]
% L is the polymer length (the number of units)
% cPolymer is the molar polymer concentration in solution
% Cmf is the molar membrane ligand concentration in solution
% fNumberOfLigands is the number of ligand types
%
% 1lig1mod: 1 ligand - 1mode of binding
% 1lig2mod: 1st ligand - modes of binding g=1 and g=2
% 1lig3mod: 1st ligand - modes of binding g=1, g=2 and g=3
%
% w(j,g1,g2) - stat weigt for a j-unit gap between g1 and g2 ligands
% w(0,:,:) - direct ligand-ligand contact
% g1-j-g2 gaps are the gaps (loops) of j units between ligands g1 and g2
% Unwrap(h,g) - the weight for a partial unwrapping of h units of g-protein
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Converted from fortran code DataModule.f90 by GJT Nov-2016

% Incorporates function DoDialogMARCKS to set up default values

function [] = SetDefaults()

global Lpolymer fNumberOfLigands rank NumberOfPoints eNumberOfChromatinStates R
global X_spacer L_spacer R_spacer
global Seq Seq1 % polymer sequence:  integer, dimension (1:10000)
global m mStart % ligand lengthes: integer, dimension (1:10) 
global mm % (the same as m(g) but real, not integer: double precision, dimension (1:10) 
global nMaxGap % max interaction lengthes : integer, dimension (1:10) 
global MAXnMaxGap MINnMaxGap % max(nMaxGap), min(nMaxGap)
global w %cooperativity parameters : double precision w(0:10000,0:10,0:10, 0:5)
global KKK % double precision(10000,10,147,5)
global Unwrap % double precision Unwrap(0:200,1:10)
global cPolymer c0 % c0(10)
global C3D C2D Cloop Alpha Beta
global productKKK
global s sigma % s(1:10000,1:5), sigma(1:5,1:5)
global iGap % states corresponding to g-j-g gaps integer, dimension (1:10000,0:10) 
global iLeftFreeEnd iRightFreeEnd % polymer end matrix state numbers
global ManualC01Range ManualC02Range ManualC03Range ManualC04Range ManualC05Range
global ManualBi CalculateMap CalculateCurves Calculate3DMap CalculateXspacer
global noGaps noLongloops
global lig1mod lig2mod lig3mod lig4mod lig5mod lig6mod

% Replacements for dialog box stuff:
global bDefaults
global Cmin Cmax
global InputFile OutputFile FastaFile SequenceFile


fNumberOfLigands = 1;
m = [2 17 17];

lig1mod = true;
lig2mod = false;
lig3mod = false;

w = zeros(7,4,4,4);

w(1,:,:,2)=[1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1];

ManualBi = false;
Cloop = 1;
Alpha = 1.75;
Beta = 0;

InputFile = 'input.dat';
OutputFile = 'output.dat';
FastaFile = 'sequence.fa';
SequenceFile = 'sequence.txt';

bDefaults = [1 0.5 0.25 0.125 0.0625 0];

nMaxGap = 0;

Lpolymer = 1000;
cPolymer = 1e-4;

c0=[0.1 0.1 0.1];

NumberOfPoints = 20;
Cmin = 1.e-12;
Cmax = 1.e-5;

noLongloops = false;

CalculateMap = false;
CalculateCurves = true;
Calculate3DMap = false;

ManualC01Range = true;
ManualC02Range = false;
ManualC03Range = false;
ManualC04Range = false;
ManualC05Range = false;




end





