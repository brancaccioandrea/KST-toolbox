function model = blim(data,fr,states,classes,options)

% BLIM   Parameter estimation of the general BLIM model
% Estimate the BLIM by maximum-likelihood via the EM algorithm.
% BLIM function allow also to put equality constraints among error
% parameters. Moreover, the EM treats missing-at-random data
%
% BLIM(DATA, FR, STATES, CLASSES, OPTIONS)
%
% estimates the parameters of the BLIM
% model for probabilistic knowledge structures by maximum likelihood,
% where:
% - data is an m-by-n binary matrix of the observed response patterns;
% - fr is a m-by-1 frequency vector of the response patterns;
% - states is a k-by-n binary matrix where each row represents a knowledge state;
% - classes is a c-by-n binary matrix specifying the equivalence classes of items. Each
%           row of classes is a class, each column an item. The matrix classes
%           must represent a partition on the set of items, therefore sum(classes)
%           must be a vector of ones and each row must contain at least a one. The
%           parameter classes can be the empty matrix []. In that case the
%           parameters of a BLIM model are estimated.
%
% BLIM(DATA,[],STATES,CLASSES)
%
% allows to pass an N-by-m binary matrix of
% ungrouped response patterns, where N is the sample size
%
% BLIM(...,OPTIONS) allows to specify the following of options:
%
%   options.xtol: tolerance for parameter change
%   options.maxiter: maximum number of iterations of EM algorithm
%   options.display: 'on' to show progress diagram, 'off' to suppress
%   options.beta0: a vector n-by-1 where each cell correspond to an item or class of item
%                  in the knowledge structure. If the value is 0 the corresponding beta parameter
%                  will not be estimated. The default is 1 on all the item
%   options.eta0:  a vector n-by-1 where each cell correspond to an item or class of item
%                  in the knowledge structure. If the value is 0 the corresponding eta parameter
%                  will not be estimated. The default is 1 on all the item
%
% Authors:  Luca Stefanutti (luca.stefanutti@unipd.it)

if nargin<4 || isempty(classes)
    classes = eye(size(states,2));
end

if nargin<5
    options.xtol=1e-5;
    options.maxiter=1000;
    options.display='off';
    options.beta0 = ones(size(classes,1),1);
    options.eta0 = ones(size(classes,1),1);
else
    if ~isfield(options,'xtol')
        options.xtol=1e-5;
    end
    if ~isfield(options,'maxiter')
        options.maxiter=1000;
    end
    if ~isfield(options,'display')
        options.display='off';
    end
    if ~isfield(options,'beta0')
        options.beta0 = ones(size(classes,1),1);
    end
    if ~isfield(options,'eta0')
        options.eta0 = ones(size(classes,1),1);
    end
end

if isempty(fr)
    [pat,fr] = patfreq(data);
else
    pat = data;
end
%
% Correct (R), wrong (W) and missing (M) parts of the response patterns
%
R = pat==1;
W = pat==0;
M = isnan(pat);
%
% Estimated probability distribution on the missing patterns
%
pm = pmiss(M,fr);

nstates=size(states,1);          % number of states
nitems=size(states,2);           % number of items
nclasses=size(classes,1);   % number of equivalence classes

sz=sum(fr); % sample size
df=2^nitems-sum(options.beta0)-sum(options.eta0)-nstates; % degrees of freedom
if df <= 0
    %    chi0=chi2inv(0.95,df); % chi-square for a 5% significance level
    %else
    chi0 = NaN;
    fprintf('\nWarning: The model has zero or negative degrees of freedom.\n');
end

chi=zeros(options.maxiter,1);
loglike=zeros(options.maxiter,1);
delta=zeros(options.maxiter,1);
exitflag=0;
%
% Initial guesses of the model's parameters
%
beta=options.beta0.*rand(nclasses,1)/10;
eta=options.eta0.*rand(nclasses,1)/10;

pk=ones(nstates,1)/nstates;
%
% main loop of the EM algorithm
%
for iter=1:options.maxiter
    %
    % expectation step
    %
    betac = classes'*beta;
    etac = classes'*eta;
    prk=rho(betac,etac,states,R,M);
    pr=prk*pk;
    prm=pm.*pr;
    pkr=(pk*(1./pr)').*prk';
    %
    % computes Log-likelihood and Pearson's chi-square statistic for
    % current iteration
    %
    loglike(iter)=-sum(fr.*log(prm));
    chi(iter)=sum(fr.*fr./(sz*prm))-sz;
    if df > 0
        pval=1-gammainc(chi(iter)/2,df/2);
    else
        pval=1;
    end
    %
    % maximization step wrt pk (the probabilities of the states)
    %
    pk_old=pk;
    pk=pkr*fr/sz;
    %
    % maximization step wrt beta and eta parameters of the item classes
    %
    beta_old = beta;
    eta_old = eta;
    for c = 1:size(classes,1)
        Kc = states(:,classes(c,:)==1);
        Rc = R(:,classes(c,:)==1);
        Wc = W(:,classes(c,:)==1);
        Mc = M(:,classes(c,:)==1);
        if options.beta0(c) > 0
            beta(c) = sum(((Kc*Wc').*pkr)*fr)/sum(((Kc*(1-Mc')).*pkr)*fr);
        end
        if options.eta0(c) > 0
            eta(c) = sum((((1-Kc)*Rc').*pkr)*fr)/sum((((1-Kc)*(1-Mc')).*pkr)*fr);
        end
    end
    %
    % Computes parameter change
    %
    theta_old = [beta_old;eta_old;pk_old];
    theta = [beta;eta;pk];
    delta(iter) = max(abs(theta_old-theta));
    %
    % plots intermediate results, if display option is on
    %
    if strcmpi(options.display,'on')==true
        figure(1);
        subplot 211
        semilogy([0,iter],[chi0,chi0],'r:');
        hold on;
        semilogy(1:iter,chi(1:iter));
        hold off;
        t=sprintf('Iteration %d. Chi2 = %6.2f. df = %d. p-value = %8.5f',iter,chi(iter),df,pval);
        title(t);
        subplot 223
        semilogy(1:iter,loglike(1:iter));
        title(sprintf('Model likelihood: %6.2f',loglike(iter)));
        subplot 224
        semilogy([0,iter],[options.xtol,options.xtol],'r:');
        hold on
        semilogy(1:iter,delta(1:iter));
        hold off
        title(sprintf('Parameter change: %f',delta(iter)));
    end
    %
    % checks if tolerance value on parameter change has been reached
    %
    if delta(iter) <= options.xtol
        %fprintf('Parameter change in EQUIBLIM was less than %f\n',options.xtol);
        exitflag=1;
        break;
    end
end
%
% Not enough iterations
%
if iter == options.maxiter
    fprintf('Maximum number of iterations exceeded in BLIM\n');
    exitflag=0;
end
%
% final computation of response pattern probabilities
%
prk=rho(classes'*beta,classes'*eta,states,R,M);
pr=prk*pk;
prm=pr.*pm;
%
% computes final Pearson's chi-square
%
model.chi=sum(fr.*fr./(sz*prm))-sz;
% if df > 0
%     model.pval=1-gammainc(model.chi/2,df/2);
% else
%     model.pval = 1;
% end
model.sz=sz;
model.df=df;
model.npar=2*nclasses+nstates-1;
model.pat=pat;
model.freq=fr;
model.expfreq=sz*prm;
model.miss=M;
model.pmiss=pm;
model.states=states;
model.classes=classes;
model.xtol=options.xtol;
model.iter=iter;
model.maxiter=options.maxiter;
model.beta=beta;
model.eta=eta;
model.pi=pk;
model.loglike=-sum(fr.*log(prm));
model.aic=2*(model.loglike+model.npar);
model.aicc=model.aic+2*model.npar*(model.npar+1)/(sz-model.npar-1);
model.bic=2*model.loglike+model.npar*log(sz);
model.exitflag=exitflag;
model.pitems=states'*pk;
if size(classes,1)~=size(classes,2)
    if sum(M,[1,2]) >0
        model.name='Costrainted-BLIM Maximum Likelihood via Expectation Maximization with Missing-at-Random';
        model.referece='de Chiusole, D., Stefanutti, L., Anselmi, P., & Robusto, E. (2018).';
    else
        model.name='Costrainted-BLIM Maximum Likelihood via Expectation Maximization';
        model.referece='de Chiusole, D., Stefanutti, L., Anselmi, P., & Robusto, E. (2018).';
    end
else
    if sum(M,[1,2]) >0
        model.name='MARBLIM Maximum Likelihood via Expectation Maximization';
        model.referece='de Chiusole, D., Stefanutti, L., Anselmi, P., & Robusto, E. (2015).';
    else
        model.name='BLIM Maximum Likelihood via Expectation Maximization';
        model.referece='Stefanutti & Robusto (2009)';
    end
end

end

function c=rho(beta,eta,states,R,M)

% Probabilistic Knowledge Structures - Local Independence Model.
% Returns the conditional probability of the patterns given
% the states.
%
% c = rho(a,b,states,patterns,M)
%
% where:
% beta    is an n-element column vector of the careless error
%          probabilities of the items;
% eta      is an n-element column vector of the lucky guess
%          probabilities of the items;
% states   is a k-by-n matrix of the states;
% R        is an m-by-n matrix of the observed response patterns;
% M
% c        is an m-by-k matrix of the conditional probabilities.
%
% Author: Luca Stefanutti (luca.stefanutti@unipd.it)

m=size(R,1);
eps=1e-9;
beta=beta+(beta<eps)*eps-(1-beta<eps)*eps;   % prevents log of zero 
eta=eta+(eta<eps)*eps-(1-eta<eps)*eps;
u=ones(m,1);
xa=u*log(beta)'; 
ya=u*log(1-beta)';
xb=u*log(eta)';
yb=u*log(1-eta)';
p=(R.*(1-M).*ya)*states';
q=(R.*(1-M).*xb)*(1-states)';
r=((1-R).*(1-M).*xa)*states';
s=((1-R).*(1-M).*yb)*(1-states)';
c=exp(p+q+r+s);

end

function [pat,freq]=patfreq(data)

% Pattern frequency function for response patterns with missing values
% (missing values encoded as NaN)

i=1; 
[m,n]=size(data);
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
        f=f+1;
        i=i+1; 
    end
    pat=[pat;p];
    freq=[freq;f];    
end

end

function pm = pmiss(M,F)

% Computation of missing pattern proportions in a data set

npat = size(M,1);
fm = zeros(npat,1);
marked = false(npat,1);
for i = 1:npat
    if ~marked(i)
        k = false(npat,1);
        for j = i:npat
            if isequal(M(i,:),M(j,:))
                k(j) = true;
            end
        end
        fm(k) = sum(F(k));
        marked = marked|k;
        if all(marked)
            break;
        end
    end
end
pm = fm/sum(F);

end
