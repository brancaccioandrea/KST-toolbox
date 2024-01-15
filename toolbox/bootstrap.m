function model0 = bootstrap(model0,nrep,options)

% BOOTSTRAP Computation of chi-square p-value via parametric bootstrap.
% 
% MODEL0 = BOOTSTRAP(MODEL0,NREP) 
% 
%  - model0 is a the outcome of a function for fitting the blim
%  - nrep the number of sample replications 
%
% BOOTSTRAP(...,OPTIONS) allows to specify the following of options:
%
%   options.xtol: tolerance for parameter change
%   options.maxiter: maximum number of iterations of EM algorithm
%   options.display: 'on' to show progress diagram, 'off' to suppress
%   options.feedback: number after how many bootstrap replications a feedback 
%                     is printed in the command window. The default is 100
%
% Authors:  Luca Stefanutti (luca.stefanutti@unipd.it)
%           Debora de Chiusole 


if nargin<3
    options.xtol=1e-4;
    options.maxiter=1000;
    options.display='off';
    options.feedback=100;
else
    if ~isfield(options,'xtol')
        options.xtol=1e-4;
    end
    if ~isfield(options,'maxiter')
        options.maxiter=1000;
    end
    if ~isfield(options,'display')
        options.display='off';
    end
    if ~isfield(options,'feedback')
        options.feedback=100;
    end
end

if ~isfield(model0,'classes')
    nitems=size(model0.states,2);
    model0.classes=eye(nitems);
end

if ~isfield(model0,'miss')
    nitems=size(model0.classes,2);
    model0.miss=zeros(model0.sz,nitems);
end

if ~isfield(model0,'LR0') && sum(model0.miss,[1,2])==0
    sat=saturated_BLIM(model0.pat,model0.freq);
    model0.saturated_loglike=sat.loglike;
    LR0=-2*(sat.loglike-model0.loglike);
elseif ~isfield(model0,'LR0')
    LR0 = nan;
else
    LR0 = model0.LR0;
end


chi0 = model0.chi;
chi = zeros(nrep,1);
pval = zeros(nrep,1);
pi=zeros(length(model0.pi),nrep);
nclasses = size(model0.classes,1);
beta = zeros(nclasses,nrep);
eta = zeros(nclasses,nrep);
iter=[];
k=0;
sz=model0.sz;

while sum(iter<model0.maxiter)<nrep
    if(length(iter)-k > nrep/5)
    fprintf('\nThe bootstrap replications reached maximum number of iterations in BLIM\n Try to reduce the tolerance value or increase the maximum number of iteration\n')
    break
    end

    [pat,freq] = sample(model0,sz);
    model = blim(pat,freq,model0.states,model0.classes,options);
    if any(isnan(pat))==false
        sat_model=saturated_BLIM(pat,freq);
    end

    if model.iter<model0.maxiter
        k=k+1;
        if rem(k,options.feedback)==0 && k/options.feedback==1
            fprintf(1,'Proportion of replication %3d%%',k/nrep*100)
        elseif rem(k,options.feedback)==0
            fprintf(1,'\b\b\b\b%3.0f%%',k/nrep*100)
        end
        beta(:,k) = model.beta;
        eta(:,k) = model.eta;
        iter(k)=model.iter;
        pi(:,k)=model.pi;
        chi(k) = model.chi;
        schi = sort(chi(1:k));
        chi_less = find(schi<=chi0);
        if isempty(chi_less)
            pc = 0;
        else
            j = max(chi_less);
            if j == k
                pc = 1;
            else
                pc = ((chi0-schi(j))/(schi(j+1)-schi(j))+j)/k;

            end
        end
        pval(k)=pc;
  
    po=(1:k)'/k;
    %pe=chi2cdf(schi,model.df);
    
    if any(isnan(pat))==false
        LR(k) = 2*(model.loglike-sat_model.loglike);
        LRs = sort(LR(1:k));
        LR_less = find(LRs<=LR0);
        if isempty(LR_less)
            pc_l = 0;
        else
            j = max(LR_less);
            if j == k
                pc_l = 1;
            else
                pc_l=((LR0-LRs(j))/(LRs(j+1)-LRs(j))+j)/k;
            end
        end
    else
        LR=nan;
        LR0=nan;
        pc_l=nan;
    end
  
    if strcmpi(options.display,'on')==true
        subplot(2,1,1);
        plot(schi,[po,pe]);
        a=axis;
        hold on;
        plot([chi0,chi0],[0,pc],'r:');
        plot([a(1),chi0],[pc,pc],'r:');
        plot(chi0,pc,'r.')
        hold off;
        title(sprintf('Replication %d. p-value = %6.4f',k,pvalue));
        xlabel('Chi-square');
        ylabel('Probability');
        subplot(2,1,2);
        plot(1:k,pval(1:k));
        drawnow;
    end
    end
end
pvalue_log = 1-pc_l;
pvalue = 1 - pc;
model0.nrep = nrep;
model0.beta_bootstrap = beta;
model0.eta_bootstrap = eta;
model0.LR_bootstrap = LR;
model0.chi_bootstrap = chi;
model0.iter_bootstrap = iter;
model0.beta_sd=std(beta,[],2,"omitnan");
model0.eta_sd=std(eta,[],2,"omitnan");
model0.beta_mean=mean(beta,2);
model0.eta_mean=mean(eta,2);
model0.pi_bootstrap=pi;
model0.pi_sd=std(pi,[],2,"omitnan");
model0.pi_mean=mean(pi,2);
model0.pvalue_chi=pvalue;
model0.pvalue_LR=pvalue_log;
model0.LR0=LR0;
end

function [pat,freq,samp] = sample(model,sz)

% SIM_EQUIBLIM   sample simulation with the EquiBLIM model

if nargin < 2
    sz = sum(model.freq);
end

nitems = size(model.states,2);

if ~isfield(model,'miss')
    model.miss=zeros(sz,nitems);
end

if ~isfield(model,'classes')
    model.classes=eye(nitems);
end

[M,fm] = patfreq(model.miss);
pm = fm/sz;

beta = model.classes'*model.beta;
eta = model.classes'*model.eta;
pc = cumsum(model.pi);
pmc = cumsum(pm);
samp = zeros(sz,nitems);
wa = model.states;
wb = 1-wa;
for i = 1:sz
    j=find(pc>=rand,1);
    m=find(pmc>=rand,1);
    a=(beta<rand(nitems,1))'.*wa(j,:);
    b=(eta>=rand(nitems,1))'.*wb(j,:);
    samp(i,:)=(a|b);
    samp(i,M(m,:)==1)=NaN;
end
[pat,freq]=patfreq(samp);
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

