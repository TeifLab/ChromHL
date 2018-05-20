% This function reads in FASTA formatted files, putting them all in
% variable Seq_DNA.

% If no arguments are passed in (nargin==0) then defaults to reading 
% all FASTA files in the directory, otherwise reads in files
% matching pattern

function [Seq_DNA] = Read_FASTA_all(varargin)


% If we have one input argument set that as pattern to read in
if (nargin==1)
    pattern = varargin{1};
else
    pattern = '*.fa*';
end

% List all files matching the pattern
FASTA_files = dir(pattern);

% If we have no files, return empty structure
if (isempty(FASTA_files))
   warning('No FASTA files found');
   Seq_DNA = [];
   return;
end

% Just select filenames in structure, no other attributes
FASTA_files = {FASTA_files.name};

% Split input pattern by directory separator
% Either "/" for UN*X or "\" for Windows
d=strsplit(pattern,filesep);

% Join d to get everything up to the last directory separator
x=length(d);

if ~isempty(d)
    d=strjoin(d(1:(x-1)),filesep);
end

% Initialise Seq_DNA object
Seq_DNA = [];

Seq_DNA = struct('file',FASTA_files,'header',cell(1,length(FASTA_files)),...
    'sequence',cell(1,length(FASTA_files)),...
    'length',cell(1,length(FASTA_files)),...
    'centre_point',cell(1,length(FASTA_files)),...
    'het_marks',cell(1,length(FASTA_files)));

% Structure of each element of Seq_DNA(i):
% file         : filename
% header       : FASTA header
% sequence     : sequence, converted to numbers [A,C,T,G,N/X,Y]=[1,2,3,4,5,6]
% length       : sequence length
% centre_point : centre point of sequence, in order to get correct number of lattice units across region
% het_marks    : predetermined heterochromatin marks, not transcription factor or repeat based

for i=1:length(FASTA_files)
    
    disp(FASTA_files{i});
    
    % Read each individual sequence: 1e6 is just greater than all sequences
    % for speed purposes as arrays are preallocated to this size in Read_FASTA
    if ~isempty(d)
        [Seq, seq_length, header] = Read_FASTA(strjoin({d,FASTA_files{i}},filesep),1e6);
    else
        [Seq, seq_length, header] = Read_FASTA(FASTA_files{i},1e6);
    end
    
    % Set up metadata for each sequence
    Seq_DNA(i).header   = header;
    Seq_DNA(i).sequence = Seq;
    Seq_DNA(i).length   = seq_length;
    Seq_DNA(i).centre_point = 20101;
    Seq_DNA(i).het_marks = zeros(1,seq_length);
    
    % As no het_marks predetermined, het_marks = 0 along all the length. Change this if
    % het_marks
    
end

end
