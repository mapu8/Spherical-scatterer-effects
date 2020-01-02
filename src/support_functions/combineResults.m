function [] = combineResults(file_prefix,N_start,N_end,outputfile)
%COMBINERESULTS Combine the calculated pressure matrixes into one matrix.

% ARGUMENTS:
% file_prefix - prefix to the filename, this should also include the path
%               to the file
% N_start and N_end - start and end values for the file suffix ((0 and 49),
%                     or (50 and 99), for example)
% outpufile - file to store the combined matrix

M = [];
for i = N_start:N_end
    filename = sprintf("%s_%d.mat",file_prefix,i);
    try
        M_temp = load(filename);
        M_temp = M_temp.M;
        M = cat(5,M,M_temp);
    catch e
        warning(sprintf("Skipped file %s",filename))
        continue
    end
end
save(outputfile, 'M');
end

