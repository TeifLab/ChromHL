tic;

job_id_string = getenv('SGE_TASK_ID');

% If it is not defined, then we are testing.
if isempty(job_id_string)
    job_id = 1;
else
    job_id = str2double(job_id_string);
end

try
    system(['mkdir -p Text_output/']);
    for kk=((job_id-1)*100):((job_id*100)-1)
        identifier=pad(num2str(kk),3,'left','0');
        Seq_DNA=Read_FASTA_all(['../individual_regions/gene_region_' identifier '*.fa']);
        if ~isempty(Seq_DNA)
            Driver;
            system(['mv output-' zz '.mat output_' sprintf('%04d',kk) '.mat']);
            system(['mv chr*.txt Text_output/']);
        end
    end
catch message
    display(['ERROR in file: ' message.stack.file])
    display(['ERROR: ' getReport(message)])
end

toc;

% Insert the matlab code you want here. Usually this will figure out
% what to do based on the job_id then start that processing.

% If we have a job ID, then exit so matlab does not sit at the command
% prompt forever.
if job_id ~= -1
    exit
end
