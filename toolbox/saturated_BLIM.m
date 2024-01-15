function model=saturated_BLIM(pat,freq)

% SATURATED_BLIM  estimates the parameters of the saturated BLIM with
% missing-at-random data. The estimation procedure is an EM algorithm
% that relies on the assumption that the observed proportion of an
% incomplete response pattern is P(R) = P(M)sum(C(R,R')P(R')) where P(M) is
% the proportion of the pattern M of missingness, P(R') is the proportion
% of the complete response pattern R' and C is the compatibility relation
% between incomplete and complete response patterns.
% 
% SATURATED_BLIM(PAT, FREQ)
%
% where:
% - pat is an m-by-n binary matrix of the observed response patterns;
% - freq is a m-by-1 frequency vector of the response patterns;
%
% Authors:  Luca Stefanutti (luca.stefanutti@unipd.it)
%           Debora de Chiusole 
maxiter = 100;
tol = 1e-4;
n=size(pat,2);

sz = sum(freq);
Rcom = expand_patterns(pat);
C = compatibility_relation(pat,Rcom);

pm = pmiss(pat,freq);
pi = ones(size(Rcom,1),1)/size(Rcom,1);
%
% Main EM loop
%
for i = 1:maxiter
    pi_old = pi;
    lambda = C'*(freq./(C*pi));
    pi = lambda.*pi/sz;
    if max(abs(pi_old-pi)) < tol
        break;
    end
end

model.pi = pi;
model.pr = C*pi;
model.pm = pm;
model.prm = model.pr.*pm;
model.loglike = -sum(freq.*log(model.prm));
model.df = 3^n-2^n;
end
function Rcom = expand_patterns(pat)
% EXPAND_PATTERNS   expands an incomplete response pattern into a matrix of
% all complete response patterns that are compatible with it.
Rcom = [];
R = pat==1;
M = isnan(pat);
nmiss = sum(M,2);
for i = 1:size(pat,1)
    C = powerset(nmiss(i));
    Ri = repmat(R(i,:),2^nmiss(i),1);
    Ri(:,M(i,:)==1) = C;
    Rcom = [Rcom;Ri];
end
Rcom = patfreq(Rcom);
end

function C = compatibility_relation(pat,Rcom)

% C(i,j) = 1 if incomplete pattern i is compatible witrh complete pattern j

R = bin2dec(num2str(pat == 1,'%d'));
M = bin2dec(num2str(isnan(pat),'%d'));
Rc = bin2dec(num2str(Rcom,'%d'));
C = zeros(size(R,1),size(Rcom,1));
for i = 1:length(R)
    U = bitor(R(i),M(i));
    C(i,:) = (bitor(R(i),Rc) == Rc) & (bitor(U,Rc) == U);
end
end

function pm = pmiss(pat,F)

% Computation of missing pattern proportions in a data set

npat = size(pat,1);
M = isnan(pat);
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
