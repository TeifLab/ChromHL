% This script is run by the SGE shell script job_submit.sh
% SGE_TASK_ID is the environment variable which holds the number of the task
% The individual FASTA files are stored in a directory ../individual_regions/gene_region_XXXX.fa
% where XXXX is a zero padded numeric ID

% Start the timer
tic;

% get the task ID from the environment variable if it exists
job_id_string = getenv('SGE_TASK_ID'); 

% If it is not defined, then we are testing.
if isempty(job_id_string)
    job_id = 1;
else % convert to a number
    job_id = str2double(job_id_string);
end

try
    % make folder to contain text output if it does not exist
    system(['mkdir -p Text_output/']);
    
    % cycle through input files
    for kk=((job_id-1)*100):((job_id*100)-1) 
        % this is kk = [0:99] + job_id*100, so one hundred runs
      
        % Set mask for identifier by padding with zeros to the left (up to three)
        identifier=pad(num2str(kk),3,'left','0');
        
        % Read the sequences matching the pattern below (in the [])
        % Change this if FASTA sequences are elsewhere
        Seq_DNA=Read_FASTA_all(['individual_regions/gene_region_' identifier '*.fa']);
        
        % If we have sequences (in the Seq_DNA structure)
        if ~isempty(Seq_DNA)
            % Pass to individual binding map calculator
            Driver; % This outputs an output file output-{zz}.mat where {zz} is a variable passed back
            
            % Change the output file in to one identified by the kk value
            system(['mv output-' zz '.mat output_' sprintf('%04d',kk) '.mat']);
            
            % Move the text output files into Text_output
            system(['mv chr*.txt Text_output/']);
        end
    end
catch message  
    % If there is an error here, retrace the stack and output all messages to stdout
    display(['ERROR in file: ' message.stack.file]) % This will show which file has the issue
    display(['ERROR: ' getReport(message)]) % This will output the exact MATLAB error found
end

toc; % Stop timer and output time elapsed - useful for timing runs and optimisation of algorithm

if job_id ~= -1 % Quit MATLAB so task ends correctly
    exit
end
