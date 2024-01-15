function info=blimit(model,filename)

% BLIMIT   Identification test of probabilistic knowledge structures
% This is a diagnostic function which is suitable for testing the
% identification of the basic local independence model for probabilistic
% knowledge structures. The identification test is based on the computation
% of the Jacobian matrix of the model prediction function. 
% Parameters of the model are: the collection of knowledge states,
% represented by a m-by-n binary matrix STATES, where m is the number of
% states and n the number of items; the careless error probabilities of the
% items, representedby a n-by-1 column vector BETA; the lucky guess
% probabilities of the items, represented by a n-by-1 column vector ETA,
% the probabilities of the knowledge states, represented by a m-by-1 column
% vector PI.
%
% BLIMIT(MODEL)
%
% prints a detailed diagnostic report on the command window, which informs
% about specific identification problems of the tested model.
%
% BLIMIT(MODEL,FILENAME)
%
% saves the report to the text file FILENAME.
%
% Author:  Luca Stefanutti (luca.stefanutti@unipd.it)

if nargin<1
    error('Not enough input arguments in call to function BLIMIT');
end
if ~isfield(model,'states')
    error('Model states are not defined');
end

states=model.states;
nstates=size(states,1);
nitems=size(states,2);

if isfield(model,'beta')
    beta=model.beta;
else
    fprintf('\nMODEL.BETA not defined: a random assignment will be done.\n');
    beta=0.49*rand(nitems,1);
end
if isfield(model,'eta')
    eta=model.eta;
else
    fprintf('\nMODEL.ETA not defined: a random assignment will be done.\n');    
    eta=0.49*rand(nitems,1);
end
if isfield(model,'pi')
    pi=model.pi;
else
    fprintf('\nMODEL.PI not defined: a random assignment will be done.\n');    
    pi=rand(nstates,1);
    pi=pi/sum(pi);
end
if isfield(model,'beta0')
    beta0=model.beta0~=0;
else
    beta0=true(nitems,1);
end
if isfield(model,'eta0')
    eta0=model.eta0~=0;
else
    eta0=true(nitems,1);
end
if nargin<2
    fid=1;
else
    fid=fopen(filename,'a');
end

% Item parameter indexes
beta_index=find(beta0);
eta_index=find(eta0);

% Tolerance for linear coefficients
eps=1e-6;

% This is the total number of parameters in the BLIM model
npar=sum(beta0)+sum(eta0)+nstates-1;

% Various arrays that are used for storing diagnostic information
% concerning the items
ta=true(1,nitems);
tb=true(1,nitems);
tab=true(2,nitems);
tap=true(1,nitems);
tbp=true(1,nitems);
tabp=true(2,nitems);

% This calls the function that computes the Jacobian matrix
jac=pksjacbas(states,beta,eta,pi);

% Decomposes the Jacobian into its main three components: 
% BETA, ETA, and PI
ja=jac(:,1:nitems);
jb=jac(:,nitems+1:2*nitems);
jp=jac(:,2*nitems+1:end);

% Removes item parameters that are fixed at zero and computes rank
ja0=ja(:,beta0);
jb0=jb(:,eta0);

jac0=[ja0,jb0,jp];
r=rank(jac0);

% Some preliminary general information
fprintf(fid,'\n\nB L I M I T\n')
fprintf(fid,'BASIC LOCAL INDEPENDENCE MODEL IDENTIFICATION ANALYSIS\n\n');
fprintf(fid,'Number of items:                 %5d\n',nitems);
fprintf(fid,'Number of knowledge states:      %5d\n\n',nstates);
fprintf(fid,'Total number of parameters:      %5d\n',npar);
fprintf(fid,'Jacobian matrix rank:            %5d\n',r);
fprintf(fid,'Null space dimension (NSD):      %5d\n\n',npar-r);

info.NItems=nitems;
info.NStates=nstates;
info.NPar=npar;
info.Rank=r;
info.NSD=npar-r;

if 2^nitems-1<=r
    fprintf(fid,'\nWARNING: The model is not quantitatively testable.\n');
    fprintf(fid,'Jacobian matrix rank (%d) is not less than the number\n',r);
    fprintf(fid,'of independent observables (%d).\n\n',2^nitems-1);
end

% Detailed identification check takes place only in case a rank deficient
% Jacobian matrix is encountered
if r<npar
    
    fprintf(fid,'Identification problems detected:\n');
    fprintf(fid,'Jacobian matrix is not full rank.\n\n');
      
    % ranks of the three submatrices
    ra=rank(ja0);
    rb=rank(jb0);
    rp=rank(jp);    
   
    % Here all pairs of submatrices are assembled and their respective
    % ranks are computed
    jab0=[ja0 jb0];
    jap0=[ja0 jp];
    jbp0=[jb0 jp];
    
    rab=rank(jab0);
    rap=rank(jap0);
    rbp=rank(jbp0);  

    % Null space dimensions
    nsda=sum(beta0)-ra;
    nsdb=sum(eta0)-rb;
    nsdp=nstates-rp-1;   
    nsdab=sum(beta0)+sum(eta0)-rab;
    nsdap=sum(beta0)+nstates-rap-1;
    nsdbp=sum(eta0)+nstates-rbp-1;
    
    % Tradeoff dimensions
    tdab=ra+rb-rab;
    tdap=ra+rp-rap;
    tdbp=rb+rp-rbp;    
    tdabp=rab+rap+rbp-ra-rb-rp-r;

    fprintf(fid,'Submatrix rank analysis table\n');
    fprintf(fid,'[BETA] = submatrix of the careless error parameters\n');
    fprintf(fid,'[ETA]  = submatrix of the lucky guess parameters\n');
    fprintf(fid,'[PI]   = submatrix of the state probabilities\n');
    fprintf(fid,'--------------------------------------------------------------\n');
    fprintf(fid,'%15s%10s%10s%10s%15s\n','SUBMATRIX','N.PAR','RANK','NSD','TRADEOFF DIM');
    fprintf(fid,'--------------------------------------------------------------\n');
    fprintf(fid,'%15s%10d%10d%10d%10d\n','[BETA]',sum(beta0),ra,nsda,nsda);
    fprintf(fid,'%15s%10d%10d%10d%10d\n','[ETA]',sum(eta0),rb,nsdb,nsdb);
    fprintf(fid,'%15s%10d%10d%10d%10d\n','[PI]',nstates-1,rp,nsdp,nsdp);
    fprintf(fid,'--------------------------------------------------------------\n');
    fprintf(fid,'%15s%10d%10d%10d%10d\n','[BETA ETA]',sum(beta0)+sum(eta0),rab,nsdab,tdab);
    fprintf(fid,'%15s%10d%10d%10d%10d\n','[BETA PI]',sum(beta0)+nstates-1,rap,nsdap,tdap);
    fprintf(fid,'%15s%10d%10d%10d%10d\n','[ETA PI]',sum(eta0)+nstates-1,rbp,nsdbp,tdbp);
    fprintf(fid,'--------------------------------------------------------------\n');
    fprintf(fid,'%15s%10d%10d%10d%10d\n','[BETA ETA PI]',npar,r,npar-r,tdabp);
    fprintf(fid,'--------------------------------------------------------------\n');    
        
    % Are there tradeoff dimensions in the [BETA] submatrix?
    if nsda>0                        
        % what the following FOR loop does is: remove from [BETA] one
        % column at the time in a stepwise fashion. Each time compute the
        % rank of the matrix and see if it remains unchanged. If so, a
        % dependent column has been found and, thus, a dependent parameter
        % has been detected
        for i=1:nitems
            index=true(nitems,1);
            index(i)=false;
            if rank(ja(:,index))==ra
                ta(i)=false;
            end
        end
        ndep=nitems-sum(ta);
        % Prints out the results of the stepwise procedure
        fprintf(fid,'\n\nItem diagnostics for [BETA] submatrix\n');
        fprintf(fid,'First-order tradeoff dimensions:   %d\n',nsda);
        fprintf(fid,'0 = parameter is not independent\n');
        fprintf(fid,'1 = parameter is independent\n');  
        fprintf(fid,'* = parameter is set to constant\n');
        fprintf(fid,'-------------------------------\n');
        fprintf(fid,'%5s%10s%15s\n','ITEMS','BETA','INDEPENDENT');
        fprintf(fid,'-------------------------------\n');
        for i=1:nitems
            fprintf(fid,'%5d%10.4f',i,beta(i));
            if beta0(i)
                fprintf(fid,'%10d\n',ta(i));
            else
                fprintf(fid,'%10s\n','*');
            end
        end
        fprintf(fid,'-------------------------------\n');
        fprintf(fid,'Lack of independence might be fixed by setting\n');
        fprintf(fid,'at zero exactly %d parameter(s) out of the %d\n',nsda,ndep);
        fprintf(fid,'parameters that are marked as not independent.\n');
    end    
    
    % Are there tradeoff dimensions in the [ETA] submatrix?
    if nsdb>0                        
        % what the following FOR loop does is: remove from [ETA] one
        % column at the time in a stepwise fashion. Each time compute the
        % rank of the matrix and see if it remains unchanged. If so, a
        % dependent column has been found and, thus, a dependent parameter
        % has been detected
        for i=1:nitems
            index=true(nitems,1);
            index(i)=false;
            if rank(jb(:,index))==rb
                tb(i)=false;
            end
        end     
        ndep=nitems-sum(tb);
        % Prints out the results of the stepwise procedure
        fprintf(fid,'\n\nItem diagnostics for [ETA] submatrix\n');
        fprintf(fid,'First-order tradeoff dimensions:   %d\n',nsdb);
        fprintf(fid,'0 = parameter is not independent\n');
        fprintf(fid,'1 = parameter is independent\n');    
        fprintf(fid,'* = parameter is set to constant\n');        
        fprintf(fid,'-------------------------------\n');
        fprintf(fid,'%5s%10s%15s\n','ITEMS','ETA','INDEPENDENT');
        fprintf(fid,'-------------------------------\n');
        for i=1:nitems
            fprintf(fid,'%5d%10.4f',i,eta(i));
            if eta0(i)
                fprintf(fid,'%10d\n',tb(i));
            else
                fprintf(fid,'%10s\n','*');
            end
        end        
        fprintf(fid,'-------------------------------\n');
        fprintf(fid,'Lack of independence might be fixed by setting\n');
        fprintf(fid,'at zero exactly %d parameter(s) out of the %d\n',nsdb,ndep);
        fprintf(fid,'parameters that are marked as not independent.\n');        
    end    
    
    % Are there tradeoff dimensions in the [BETA ETA] submatrix?
    if tdab>0                        
        % what the following FOR loop does is: remove from [BETA ETA] one
        % column at the time in a stepwise fashion. Each time compute the
        % rank of the matrix and see if it remains unchanged. If so, a
        % dependent column has been found and, thus, a dependent parameter
        % has been detected
        for i=1:nitems  
            index=true(nitems,1);            
            index(i)=false;
            if rank([ja(:,index) jb])==rab
                tab(1,i)=false;
            end
        end
        ndep=nitems-sum(tab(1,:));
        for i=1:nitems  
            index=true(nitems,1);            
            index(i)=false;
            if rank([ja jb(:,index)])==rab
                tab(2,i)=false;
            end
        end
        ndep=ndep+nitems-sum(tab(2,:));        
        % Prints out the results of the stepwise procedure
        fprintf(fid,'\n\nItem diagnostics for [BETA ETA] submatrix\n');
        fprintf(fid,'Second-order tradeoff dimensions:   %d\n',tdab);
        fprintf(fid,'0 = parameter is not independent\n');
        fprintf(fid,'1 = parameter is independent\n');       
        fprintf(fid,'* = parameter is set to constant\n');                
        fprintf(fid,'--------------------------------------------------------\n');
        fprintf(fid,'%5s%10s%15s','ITEMS','BETA','INDEPENDENT');
        fprintf(fid,'%10s%15s\n','ETA','INDEPENDENT');
        fprintf(fid,'--------------------------------------------------------\n');
        for i=1:nitems
            fprintf(fid,'%5d%10.4f',i,beta(i));
            if beta0(i)
                fprintf(fid,'%10d',tab(1,i));
            else
                fprintf(fid,'%10s','*');
            end
            fprintf(fid,'%15.4f',eta(i));
            if eta0(i)
                fprintf(fid,'%10d\n',tab(2,i));
            else
                fprintf(fid,'%10s\n','*');
            end            
        end      
        fprintf(fid,'--------------------------------------------------------\n');
        fprintf(fid,'Lack of independence might be fixed by setting\n');
        fprintf(fid,'at zero exactly %d parameter(s) out of the %d\n',tdab,ndep);
        fprintf(fid,'parameters that are marked as not independent.\n');        
    end
    
    % Is the [BETA PI] submatrix full rank? If not then there are tradeoffs
    % between some BETA and some PI parameters
    beta0_new=beta0;
    if tdap>0          
        nsb=null([jp ja0],'r');
        fprintf(fid,'\n\nItem diagnostics for [BETA PI] submatrix\n');
        fprintf(fid,'Second-order tradeoff dimensions:   %d\n',tdap);               
        printline(fid,16+10*size(nsb,2)+1);
        fprintf(fid,'%6s%10s','PARAMS','VALUES');
        for j=1:size(nsb,2)
            hdr=sprintf('DIM #%d',j);
            fprintf(fid,'%10s',hdr);
        end
        fprintf(fid,'\n');
        printline(fid,16+10*size(nsb,2)+1);
        for i=1:length(beta_index)
            par=sprintf('%4s%d','BETA',beta_index(i));
            fprintf(fid,'%6s%10.4f',par,beta(beta_index(i)));
            for j=1:size(nsb,2)
                coeff=nsb(nstates+i-1,j);
                if abs(coeff)>eps
                    fprintf(fid,'%10.4f',coeff);                      
                else
                    fprintf(fid,'%10s',' ');
                end
                if coeff==1
                    beta0_new(beta_index(i))=false;
                end
            end
            fprintf(fid,'\n');
        end
        printline(fid,16+10*size(nsb,2)+1);
        tap=nsb;        
    end
    
    % Is the [ETA PI] submatrix full rank? If not then there are tradeoffs
    % between some ETA and some PI parameters
    eta0_new=eta0;
    if tdbp>0         
        nsb=null([jp jb0],'r');
        fprintf(fid,'\n\nItem diagnostics for [ETA PI] submatrix\n');
        fprintf(fid,'Second-order tradeoff dimensions:   %d\n',tdbp);               
        printline(fid,16+10*size(nsb,2)+1);
        fprintf(fid,'%6s%10s','PARAMS','VALUES');
        for j=1:size(nsb,2)
            hdr=sprintf('DIM #%d',j);
            fprintf(fid,'%10s',hdr);
        end
        fprintf(fid,'\n');
        printline(fid,15+10*size(nsb,2)+1);
        for i=1:length(eta_index)
            par=sprintf('%4s%d','ETA',eta_index(i));
            fprintf(fid,'%6s%10.4f',par,eta(eta_index(i)));
            for j=1:size(nsb,2)
                coeff=nsb(nstates+i-1,j);
                if abs(coeff)>eps
                    fprintf(fid,'%10.4f',coeff);                      
                else
                    fprintf(fid,'%10s',' ');
                end
                if coeff==1
                    eta0_new(eta_index(i))=false;                    
                end
            end
            fprintf(fid,'\n');
        end
        printline(fid,16+10*size(nsb,2)+1);        
        tbp=nsb;
    end

    % So far only tradeoffs between pairs of parameters have been checked.
    % Here we check whether there are tradeoffs involving all three types
    % of parameters. This happens when all three submatrices [BETA ETA],
    % [BETA PI] and [ETA PI] are full rank while [BETA ETA PI] is not.
    if tdabp>0
        for i=1:nitems
            index=true(nitems,1);
            index(i)=false;
            if rank([ja(:,index) jb jp])==r
                tabp(1,i)=0;
            end
        end
        for i=1:nitems
            index=true(nitems,1);
            index(i)=false;
            if rank([ja jb(:,index) jp])==r
                tabp(2,i)=0;
            end
        end     
        ja0_new=ja(:,beta0_new);
        jb0_new=jb(:,eta0_new);
        beta_index_new=find(beta0_new);
        eta_index_new=find(eta0_new);
        nsb=null([jp ja0_new jb0_new],'r');        
        fprintf(fid,'\n\nItem diagnostics for [BETA ETA PI] submatrix\n');
        fprintf(fid,'Third-order tradeoff dimensions:   %d\n',tdabp);               
        printline(fid,16+10*size(nsb,2)+1);
        fprintf(fid,'%6s%10s','PARAMS','VALUES');
        for j=1:size(nsb,2)
            hdr=sprintf('DIM #%d',j);
            fprintf(fid,'%10s',hdr);
        end
        fprintf(fid,'\n');
        printline(fid,16+10*size(nsb,2)+1);
        for i=1:length(beta_index_new)
            par=sprintf('%4s%d','BETA',beta_index_new(i));
            fprintf(fid,'%6s%10.4f',par,beta(beta_index_new(i)));
            for j=1:size(nsb,2)
                coeff=nsb(nstates+i-1,j);
                if abs(coeff)>eps
                    fprintf(fid,'%10.4f',coeff);               
                else
                    fprintf(fid,'%10s',' ');
                end               
            end
            fprintf(fid,'\n');
        end
        printline(fid,16+10*size(nsb,2)+1);        
        for i=1:length(eta_index_new)
            par=sprintf('%4s%d','ETA',eta_index_new(i));
            fprintf(fid,'%6s%10.4f',par,eta(eta_index_new(i)));
            for j=1:size(nsb,2)
                coeff=nsb(nstates+length(beta_index_new)+i-1,j);
                if abs(coeff)>eps
                    fprintf(fid,'%10.4f',coeff);               
                else
                    fprintf(fid,'%10s',' ');
                end              
            end
            fprintf(fid,'\n');
        end
        printline(fid,16+10*size(nsb,2)+1);   
        tabp=nsb;
    end
    
    info.RankBeta=ra;
    info.RankEta=rb;
    info.RankPi=rp;
    info.RankBetaEta=rab;
    info.RankBetaPi=rap;
    info.RankEtaPi=rbp; 

else
    fprintf(fid,'No identification problems detected.\n\n');
end

info.DiagBetaEta=tab;
info.DiagBetaPi=tap;
info.DiagEtaPi=tbp;
info.DiagBetaEtaPi=tabp;
info.Jacobian=jac;

if fid~=1
    fclose(fid);
    fprintf('\n\nResults saved in %s\n\n',filename);
end


function [jac,basis]=pksjacbas(states,beta,eta,pi)

% PKSJACBAS   Computation of the Jacobian matrix of the log-likelihood
% function of the basic local independence model for probabilistic
% knowledge structures. 
%
% JAC=PKSJACBAS(STATES,ALPHA,BETA,PI)
%
% For m the number of response patterns, and n the number of free
% parameters, if m > n, PKSJACBAS does not compute the full Jacobian
% matrix. The returned matrix contains at most n rows, representing a
% minimal spanning set (basis) for the vector space of the full Jacobian.
% If the Jacobian is not full rank, typically this minimal spanning set
% will have less than n vectors.
% Input arguments: a binary matrix STATES of the knowledge states; item
% parameter vectors ALPHA and BETA; a vector of the states probabilities
% PI. Arguments ALPHA, BETA and PI are optional. 
% Output arguments: the n by n Jacobian matrix JAC; a binary matrix BASIS
% whose rows are the response patterns corresponding to the rows of the
% matrix JAC.

if nargin<1
    error('Not enough input arguments');
end

[nstates,nitems]=size(states);

% this is the total number of free parameters in the model
npar=2*nitems+nstates-1;    

if nargin<2 || isempty(beta)
    beta=0.5*rand(nitems,1);
else
    if size(beta,1)~=nitems || size(beta,2)~=1
        error('Wrong size of beta vector');
    end
end

if nargin<3 || isempty(eta)
    eta=0.5*rand(nitems,1);
else
    if size(eta,1)~=nitems || size(eta,2)~=1
        error('Wrong size of eta vector');
    end
end

if nargin<4 || isempty(pi)
    pi=rand(nstates,1);
    pi=pi/sum(pi);
else
    if size(pi,1)~=nstates || size(pi,2)~=1
        error('Wrong size of pi vector');
    end
end

% initial basis is the empty set, initial Jacobian is the empty matrix,
% initial rank is zero
basis=[];
jac=[];
r=0;

% main loop executes for at most 2^nitems-1 times
for d=0:2^nitems-1
    % convert integer scalar d into binary response pattern p
    x=d;
    k=nitems;
    p=zeros(1,nitems);
    while k>0 && x>0
        p(k)=rem(x,2);
        x=floor(x/2);
        k=k-1;
    end

    % compute Jacobian row corresponding to p and
    % add it, as a new line, at the bottom of matrix jac. 
    % If rank(jac) increases then retain p, otherwise discard it.
    basnew=[basis;p];
    jacrow=pksjac(p,states,beta,eta,pi);
    jacnew=[jac;jacrow];
    rnew=rank(jacnew);
    if rnew>r
        basis=basnew;
        jac=jacnew;
        r=rnew;
        fprintf('\n%d/%d:\t',r,npar);        
        fprintf('%d',p);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fprintf('\t %3.3f',jacrow);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    if size(basis,1)==npar % ...we are done
        break;
    end
end
fprintf('\n\n');


function jac=pksjac(patterns,states,a,b,pi)

% Computation of the Jacobian matrix of the BLIM model

eps=1e-9;
m=size(patterns,1);
a=a+(a<eps)*eps-(1-a<eps)*eps;  
b=b+(b<eps)*eps-(1-b<eps)*eps;
u=ones(m,1);
xa=u*log(a)'; 
ya=u*log(1-a)';
xb=u*log(b)';
yb=u*log(1-b)';
p=(patterns.*ya)*states'+(patterns.*xb)*(1-states)';
p=p+((1-patterns).*xa)*states'+((1-patterns).*yb)*(1-states)';
c=exp(p);
xa=c*diag(pi)*states;
xb=c*diag(pi)*(1-states);
da=(xa.*(1-patterns))*diag(1./a)-(xa.*patterns)*diag(1./(1-a));
db=(xb.*patterns)*diag(1./b)-(xb.*(1-patterns))*diag(1./(1-b));
dp=c-repmat(c(:,end),1,size(c,2));
dp=dp(:,1:size(c,2)-1);
jac=[da,db,dp];


function printline(fid,length)

for i=1:length
    fprintf(fid,'-');
end
fprintf(fid,'\n');