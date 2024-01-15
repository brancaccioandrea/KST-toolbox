function [kspace,cspace,itemindex]=skillmap(map,cspace)
% SKILLMAP   Derivation of knolwedge strucuture from a skill map
% Return the knowledge structure obtained from a skill-map
%
% KSPACE=SKILLMAP(MAP)
%
% [KSPACE,CSPACE,ITEMINDEX]=SKILLMAP(MAP)
%
% where:
% - map consists of a item indicator, in the forms of a numeric vectore and a
%       items by skill matrix where a skills j is require from item j.
% - cspace the competence space. The default containting 2^|skills| comepetence states
%
% Author:  Luca Stefanutti (luca.stefanutti@unipd.it)
%
map=sortrows(map,1);
index=zeros(size(map,1),1);

k=1;
index(1)=map(1,1);
for i=2:size(map,1)
    if map(i,1)>map(i-1,1)
        k=k+1;        
    end
    index(i)=k;
end
nitems=max(index);
nskills=size(map,2)-1;

if nargin<2
    cspace=powerset(nskills);
end
nstates=size(cspace,1);
kspace=zeros(nstates,nitems);

for k=1:nstates
    for i=1:size(map,1)
        x=map(i,2:nskills+1);
        if (x&cspace(k,:))==x
            kspace(k,index(i,1))=1;
        end
    end
end

itemindex=zeros(nitems,1);
for i=1:nitems
    j=find(index(:,1)==i,1,'first');
    itemindex(i)=map(j,1);
end
    

    