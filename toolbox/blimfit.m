function best_model=blimfit(data,states,nrep_gm,nrep_bootstrap,classes,options)
% BLIMFIT Parameter estimation of the general BLIM model with bootstrap
% p-value.Estimate the BLIM by maximum-likelihood via the EM algorithm.
% BLIMFIT function allow also to put equality constraints among error
% parameters. Moreover, the EM treats missing-at-random data
%
% BEST_MODEL=BLIMFIT(DATA, STATES,NREP,NREP_BOOTSTRAP,CLASSES, OPTIONS)
%
% estimates the parameters of the BLIM
% model for probabilistic knowledge structures by maximum likelihood,
% where:
% - data is an N-by-m binary matrix of ungrouped response patterns, 
%        where N is the sample size
% - states is a k-by-n binary matrix where each row represents a knowledge state;
% - nrep_gm is the number of times in which the data are fitted starting from different points
%           of the parameter space. Thethe replication obtaining the highest model likelihood is selected.
% - nrep_bootstrap is the number of sample replications involved to compute the bootstrap p-value. 
% - classes is a c-by-n binary matrix specifying the equivalence classes of items. Each
%           row of classes is a class, each column an item. The matrix classes
%           must represent a partition on the set of items, therefore sum(classes)
%           must be a vector of ones and each row must contain at least a one. The
%           parameter classes can be the empty matrix []. In that case the
%           parameters of a BLIM model are estimated.
%
% BLIMFIT(...,OPTIONS) allows to specify the following of options:
%
%   options.xtol: tolerance for parameter change
%   options.maxiter: maximum number of iterations of EM algorithm
%   options.display: 'on' to show progress diagram, 'off' to suppress
%   options.nitems: defines the maximum number of items for which the LR 
%                   statistic is computed. The default is 15 
%
% Authors:  Andrea Brancaccio (andrea.brancaccio@unipd.it)

if nargin<4 
    nrep_bootstrap = 0;
end

if nargin<5 || isempty(classes)
    classes = eye(size(states,2));
end

if nargin<6
    options.tol = 1e-4;
    options.maxiter = 1000;
    options.display = 'off';
    options.nitems = 15;
else
    if ~isfield(options,'tol')
        options.tol = 1e-4;
    end
    if ~isfield(options,'maxiter')
        options.maxiter = 1000;
    end
    if ~isfield(options,'display')
        options.display = 'off';
    end 
    if ~isfield(options,'nitems')
        options.nitems = '15';
    end 
end

[pat,fi]=patfreq(data);
model=cell(nrep_gm,1);
loglike=zeros(nrep_gm,1);

% Fitting the model 'nrep' times to avoid local minimum
fprintf('\nLocal minima checking...\n')
for i=1:nrep_gm
    model{i,1}=blim(pat,fi,states,classes,options);
    loglike(i)=model{i,1}.loglike;
end

% Selecting the model having the smallest loglikelihood
[~,best_id]=min(loglike);
best_model=model{best_id,1};

% Computing the likelihood ratio test (LR) with the saturated BLIM
if options.nitems>=size(states,2)
    sat=saturated_BLIM(pat,fi);
    best_model.saturated_loglike=sat.loglike;
    best_model.LR0=-2*(sat.loglike-best_model.loglike);
else
    fprintf('\nThe likelihood ratio (LR) test will not be compute\nIncrease the options.nitems argument to compute the LR test\n\n')
end

% Bootstrapped p-value computation
if nrep_bootstrap > 0
    [best_model]=bootstrap(best_model,nrep_bootstrap,options);
end
printmodel(best_model)
end



function [pat,freq]=patfreq(data,fr)

% Pattern frequency function for response patterns with missing values
% (missing values encoded as NaN)

i=1;
[m,n]=size(data);

if nargin < 2
    fr=ones(m,1);
end

pat=[];
freq=[];
data1=data;
data1(isnan(data1))=2;
b=3.^(n-1:-1:0);
d=data1*b';
[sd,j]=sort(d+1);
while i<=length(sd)
    f=0;
    v=sd(i);
    p=data(j(i),:);
    while i<=length(sd)&&sd(i)==v
        f=f+fr(j(i));
        i=i+1;
    end
    pat=[pat;p];
    freq=[freq;f];
end
end
