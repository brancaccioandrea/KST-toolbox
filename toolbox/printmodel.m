function printmodel(model,filename)
% PRINTMODEL Print in console or on file the summary results for the BLIM 
% Return the estimates and the statistics for the BLIM and its extension
%
% PRINTMODEL(MODEL,FILENAME,MODE)
%
% where:
% - model is the outcome of an application of the BLIM
% - filename is the name and paths of the file in which you  want to save
%   your output. For instance '\PATH\output.txt'. If filename is empty then
%   the outcome is printed in console
%
% Authors:  Andrea Brancaccio (andrea.bracaccio@unipd.it)
%           Luca Stefanutti 
%           Debora de Chiusole 

if nargin>1
    fid=fopen(filename,mode);
else
    fid=1;
end
fprintf(fid,'\n*** Basic Local Independence Model ***\n');
if isfield(model,'name')
    fprintf(fid,'Title: %s\n\n',model.name);
else
    fprintf(fid,'\n');
end
if isfield(model,'referece')
    fprintf(fid,'Reference: %s\n\n',model.referece);
else
    fprintf(fid,'\n');
end
fprintf(fid,'Knolwedge structure\n');
fprintf(fid,'\tNumber of items:               %d\n',size(model.states,2));
fprintf(fid,'\tNumber of knowledge states:    %d\n',size(model.states,1));

fprintf(fid,'\nSample\n');
fprintf(fid,'\tNumber of observations:         %d\n',model.sz);
fprintf(fid,'\tNumber of response patterns:    %d\n',size(model.pat,1));
fprintf(fid,'\nAlgorithm: EM\n');
fprintf(fid,'\tConvergence criterion (parameters):%f\n',model.xtol);
fprintf(fid,'\tNumber of maximum iterations:   %d\n',model.maxiter);
fprintf(fid,'\tActual iteration:              %d\n',model.iter);
if model.exitflag==1
    fprintf(fid,'\nThe algorithm terminates after reaching the tolerance value\n\n');
    fprintf(fid,'Fit of the model:\n');
    fprintf(fid,'\tChi-Squared:                  %6.2f\n',model.chi);
    fprintf(fid,'\tDegrees of freedom:             %d\n',model.df);
    if isfield(model,'nrep')
        fprintf(fid,'\tBootstrapped p-value           %f\n',model.pvalue_chi);
        fprintf(fid,'\tBootstrap replications         %i\n',model.nrep);
     if ~isnan(model.LR0)
        fprintf(fid,'\nComparison with the Saturated Model:\n');
        fprintf(fid,'\tLoglikelihood of the model:    %6.2f\n',-model.loglike);
        fprintf(fid,'\tLikelihood ratio:             %6.2f\n',model.LR0);
        fprintf(fid,'\tBootstrapped p-value           %f\n',model.pvalue_LR);
        fprintf(fid,'\tBootstrap replications         %i\n',model.nrep);
     end
    end
    fprintf(fid,'\n');
    fprintf(fid,'\nParameter estimates\n\n');
    if isfield(model,'classes')
        [class,item]=find(model.classes);
    fprintf(fid,'%5s%5s%13s%13s\n','class','item','beta','eta');
        for i=1:size(model.classes,2)
            if model.eta(class(i))+model.beta(class(i))>=1
                fprintf(fid,'%5d%5d%13f%13f  **Warnings beta+eta>1\n',class(i),item(i),model.beta(class(i)),model.eta(class(i)));
            elseif model.beta(class(i))>=.5
                fprintf(fid,'%5d%5d%13f%13f  **Warnings beta>.5\n',class(i),item(i),model.beta(class(i)),model.eta(class(i)));
            elseif  model.eta(class(i))>=.5
                fprintf(fid,'%5d%5d%13f%13f  **Warnings eta>.5\n',class(i),item(i),model.beta(class(i)),model.eta(class(i)));
            else
                fprintf(fid,'%5d%5d%13f%13f\n',class(i),item(i),model.beta(class(i)),model.eta(class(i)));
            end
        end
    end
else
    fprintf(fid,'\nThe algorithm terminates after reaching the maximum number of iteration\n\n');
end


if nargin==2
    fclose(fid);
end