function [Seq_DNA] = Read_FASTA_all(varargin)


if (nargin==1)
    pattern = varargin{1};
else
    pattern = '*.fa*';
end

FASTA_files = dir(pattern);

if (isempty(FASTA_files))
   warning('No FASTA files found');
   Seq_DNA = [];
   return;
end

FASTA_files = {FASTA_files.name};

d=strsplit(pattern,filesep);


x=length(d);

if ~isempty(d)
    d=strjoin(d(1:(x-1)),filesep);
end

%disp(d)

Seq_DNA = [];

Seq_DNA = struct('file',FASTA_files,'header',cell(1,length(FASTA_files)),...
    'sequence',cell(1,length(FASTA_files)),...
    'length',cell(1,length(FASTA_files)),...
    'centre_point',cell(1,length(FASTA_files)),...
    'het_marks',cell(1,length(FASTA_files)));

    for i=1:length(FASTA_files)
    
    disp(FASTA_files{i});
    
    if ~isempty(d)
        [Seq, seq_length, header] = Read_FASTA(strjoin({d,FASTA_files{i}},filesep),1e6);
    else
        [Seq, seq_length, header] = Read_FASTA(FASTA_files{i},1e6);
    end
    
    Seq_DNA(i).header   = header;
    Seq_DNA(i).sequence = Seq;
    Seq_DNA(i).length   = seq_length;
    Seq_DNA(i).centre_point = 20101;
    Seq_DNA(i).het_marks = zeros(1,seq_length);
    %    Seq_DNA(i).het_marks(21) = 1;
    
end


end
