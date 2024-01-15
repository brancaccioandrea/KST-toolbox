function [P]=powerset(n)
% POWERSET  Computes the power set of a set
% Return the powerset  of cardinalities 2^n
%
% P = POWERSET(n)
%
% where:
% - n in the number of items in the domain of knolwedge
%
% Authors:  Luca Stefanutti (luca.stefanutti@unipd.it)
%           Debora de Chiusole

k=2^n;
P=zeros(k, n);
for i=1:k
	P(i,:)=dec2bin(i-1,n);
end

function [b]=dec2bin(d,n)

dec=d;
k=n;
b=zeros(1,n);
while(k > 0 && dec > 0)
	b(k)=rem(dec,2);
	dec=floor(dec/2);
	k=k-1;
end
