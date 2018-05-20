% This file calculates the binding map (only) based on
% input sequence
%
% It needs, at least, the DNA_sequence, converted from ACGT to 1234
% and pre-determined heterochromatin marks
%
% If a third argument is given, this is the fixed-width text file for output
%
function [profile] = CalculateMap(varargin)

global Lpolymer fNumberOfLigands rank NumberOfPoints eNumberOfChromatinStates R
global X_spacer L_spacer R_spacer
global Seq Seq1 % polymer sequence:  integer, dimension (1:10000)
global m mStart % ligand lengthes: integer, dimension (1:10)
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

% Check for numbers of arguments in 
if (nargin==3)
    DNA_seq = varargin{1};
    heterochromatin_marks = varargin{2};
    fileout = varargin{3};
elseif (nargin==2)
    DNA_seq = varargin{1};
    heterochromatin_marks = varargin{2};
    fileout = '';
end

% Set default values of parameters (if necessary)
SetDefaults;

% Set fixed-width formats for output
% Headers:
format3columns = '%4s %19s %19s\n';
format4columns = '%4s %19s %19s %19s\n';
format5columns = '%4s %19s %19s %19s %19s\n';
format6columns = '%4s %19s %19s %19s %19s %19s\n';
format7columns = '%4s %19s %19s %19s %19s %19s %19s\n';
format8columns = '%4s %19s %19s %19s %19s %19s %19s %19s\n';
format9columns = '%4s %19s %19s %19s %19s %19s %19s %19s %19s\n';
format10columns = '%4s %19s %19s %19s %19s %19s %19s %19s %19s %19s\n';
format11columns = '%4s %19s %19s %19s %19s %19s %19s %19s %19s %19s %19s\n';
format12columns = '%4s %19s %19s %19s %19s %19s %19s %19s %19s %19s %19s %19s\n';
format13columns = '%4s %19s %19s %19s %19s %19s %19s %19s %19s %19s %19s %19s %19s\n';
format14columns = '%4s %19s %19s %19s %19s %19s %19s %19s %19s %19s %19s %19s %19s %19s\n';

% Output_data formats:
format_I2F = '%4d%20.10E%20.10E\n';
format_I3F = '%4d%20.10E%20.10E%20.10E\n';
format_I4F = '%4d%20.10E%20.10E%20.10E%20.10E\n';
format_I5F = '%4d%20.10E%20.10E%20.10E%20.10E%20.10E\n';
format_I6F = '%4d%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E\n';
format_I7F = '%4d%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E\n';
format_I8F = '%4d%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E\n';
format_I9F = '%4d%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E\n';
format_I10F = '%4d%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E\n';
format_I11F = '%4d%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E\n';
format_I12F = '%4d%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E\n';
format_I13F = '%4d%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E%20.10E\n';

% nMaxGap holds the size of the maximum gap between proteins
nMaxGap = nMaxGap*ones(fNumberOfLigands,1);

% If no gaps between ligands is possible, set the noGaps flag
if(nMaxGap(1)<=1)
    noGaps=true;
else
    noGaps=false;
end

% Fill in the parameters for the simulation
ParametersInitMicrodomain(DNA_seq,heterochromatin_marks); % fill in K parameters, s, sigma
ConstantsInitMicrodomain; % fill in remaining Unwrap values &c


% Do the binding map computation
fid8 = -1;

% if a text file is given, open it for writing
if ~isempty(fileout)
    fid8=fopen(fileout,'w');
end

% Text file in output has header of form:
%
% header is Site, Ci{1...f}, theta{1...e}
%
% then the file contains
% Site     site number, numbered from left edge (1)
% Ci[j]    is the binding probability of ligand j
% theta[j] is the probability of heterochromatin state j

% This is total number of output columns
zTotalOutputs=fNumberOfLigands+eNumberOfChromatinStates;

% Calculate map
[cMap,tetaMap]=MapOfBindingCalc;

formats={format3columns,format4columns,...
    format5columns,format6columns,format7columns,...
    format8columns,format9columns,format10columns,...
    format11columns,format12columns,format13columns,...
    format14columns};

header=cell(1+zTotalOutputs,1);

header{1}='Site';
for i=1:fNumberOfLigands
    header{i+1}=strcat('Ci',num2str(i));
end
for j=1:eNumberOfChromatinStates
    header{j+1+fNumberOfLigands}=strcat('theta',num2str(j));
end

profile = [cMap(2:fNumberOfLigands+1,:); tetaMap(1:eNumberOfChromatinStates,:)];

% Write header to file
if (and(~isempty(fileout),fid8>0))
    if zTotalOutputs<=13
        fprintf(fid8, formats{zTotalOutputs-1}, header{:});
    else
        format = [formats{12} repmat(' %19s',zTotalOutputs-13,1)];
        fprintf(fid8, format, header{:});
    end
end

formats={format_I2F,format_I3F,format_I4F,format_I5F,format_I6F,...
    format_I7F,format_I8F,format_I9F,format_I10F,format_I11F,...
    format_I12F,format_I13F};

% Output binding map line-by-line in fixed format
for TestSiteNumber=1:Lpolymer

    if (and(~isempty(fileout),fid8>0))
        if(zTotalOutputs<=13)
            fprintf(fid8,formats{zTotalOutputs-1}, TestSiteNumber,...
                cMap(2:fNumberOfLigands+1,TestSiteNumber),...
                tetaMap(1:eNumberOfChromatinStates,TestSiteNumber));
        else
            format = [formats{12} repmat('%20.10E',zTotalOutputs-13,1)];
            fprintf(fid8,format, TestSiteNumber,...
                cMap(2:fNumberOfLigands+1,TestSiteNumber),...
                tetaMap(1:eNumberOfChromatinStates,TestSiteNumber));
        end
    end
    
end

% Close file if opened before
if (fid8>0) 
    fclose(fid8);
end

% return binding map in variable profile
end

