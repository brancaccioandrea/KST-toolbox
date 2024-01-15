function states=base2struct(basis)
% BASE2STRUCT Building a knowledge space from its base. 
% Implementation of the algorithm described in Doignon & Falmagne, 1999, p.32.
% Return the knowledge strucuture
%
% STATES=BASE2STRUCT(BASIS)
%
% where:
% - basis is a atoms by items dicothoums matrix
% - states is a knowledge states by items matrix where the elements in a
%         (i,j) cell is equal to 1 if the items j belong to knowledge state i, otherwise is 0.
%
% Authors:  Luca Stefanutti (luca.stefanutti@unipd.it)
%           Debora de Chiusole 
%
natoms=size(basis,1);
nitems=size(basis,2);
emptyset=zeros(1,nitems);
states=emptyset;
for i=1:natoms
    h=[];
    for j=1:size(states,1)
        if ~subset(basis(i,:),states(j,:))
            flag=true;
            u=basis(i,:)|states(j,:);
            for k=1:i-1                
                if subset(basis(k,:),u)
                    if ~subset(basis(k,:),states(j,:))
                        flag=false;
                        break;
                    end
                end
            end
            if flag==true
                h=[h;u];
            end
        end
    end
    states=[states;h];
end

function retval=subset(x,y)

% return true if x is a subset of y

retval=isequal(x&y,x);




