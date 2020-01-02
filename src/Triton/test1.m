function [] = test1(outputfile,N_diffuse_reps)
%TEST1 Summary of this function goes here
%   Detailed explanation goes here
try
    outputfile
    N_diffuse_reps
catch e
    fprintf(1,'The identifier was:\n%s\n',e.identifier);
    fprintf(1,'There was an error! The message was:\n%s',e.message);
    exit(0)
end
end

