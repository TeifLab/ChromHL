% This function reads in an individual FASTA file, given in
% FastaFile, with expected length Length

function [Seq,seq_length,header] = Read_FASTA(FastaFile,Length)

Seq = zeros(1,Length);

% Try opening the FASTA file
fid4 = fopen(FastaFile,'r');

seq_length = 0;

if fid4>0 % File exists and can be opened
    
    % Get first line of file
    header = fgets(fid4);
    
    % If doesn't start with ">", then not a FASTA file
    if header(1)~='>'
        warning(['''' FastaFile ''' is not a FASTA file']);
    else
        % Strip off the first character, and write in header variable
        header=header(2:end-1);
        
        % If FASTA file has a CR at the end (Windows), strip that off
        if uint8(header(end))==13
            % (with CRLF rather than just LF)
            header=header(1:end-1);
        end
    
        disp(['Sequence for ' header ' found']);
        i = 1;
        
        % Cycle along sequence, reading sequence as we go
        while i<=Length
            % Scan character by character
            nextsymbol = fscanf(fid4,'%c',1);
            
            % If at end of file, break loop 
            if ~ischar(nextsymbol) || isempty(nextsymbol)
                warning(['Only ' num2str(i-1) ' symbols found in sequence']);
                break;
            end
            
            % If character is LF, CR or space move to next character
            % (FASTA file can be fixed-width format, so this just skips to next symbol)
            if uint8(nextsymbol)==10 || uint8(nextsymbol)==13  || uint8(nextsymbol)==32
                continue;
            else % we have a (suspected) nucleotide
                switch lower(nextsymbol) % into lower-case
                    case 'a'
                        Seq(i)=1;
                    case 'c'
                        Seq(i)=2;
                    case 'g'
                        Seq(i)=3;
                    case {'t','u'} % Could be thyamine (T) or uracil (U) if RNA seq given
                        Seq(i)=4;
                    case {'x','n'} % If (N) or (X)
                        Seq(i)=5;
                    case 'y'
                        Seq(i)=6;
                    otherwise % We don't have a nucleotide, possibly amino acid symbol?
                        warning(['Unknown sequence symbol: ' nextsymbol]);
                        % keep going, could just be a gap
                end
            end
            % Move onto next symbol
            i=i+1;
        end
  
        % Now seq_length contains number of symbols in sequence+1
        seq_length = i-1;
    end
else
    seq_length = 0;
    disp(['''' FastaFile ''' is not readable']);
end

% Trim off trailing 0s in Seq variable, as not pertinent
Seq = Seq(1:seq_length);

fclose(fid4);

end
