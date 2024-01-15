function hGraphPlot = plot_knowledge_structure(states,labelstyle,labels)
% PLOT_KNOWLEDGE_STRUCTURE Plots the Hasse diagram of a knowledge structure
% Return a directed graph object 
%
% HGRAPHPLOT=PLOT_KNOWLEDGE_STRUCTURE(states,labelstyle)
%
% where:
% - states is a knowledge states by items matrix where the elements in a
%          (i,j) cell is equal to 1 if the items j belong to knowledge state i,
%          otherwise is 0.
% -labelstyle is a string determined the type of label associate with each
%          knolwedge state in the graph. Labelstyle variable may be: 
%          numbers
%          vectors
%          literals
%          cardinality
%          custom
% - label is a vector of string of the length of the number of item and it
%         is used for the 'custom' labelstyle 
%
%
% Authors:  Luca Stefanutti (luca.stefanutti@unipd.it)
%           Debora de Chiusole 
% 

if nargin < 2
    labelstyle = 'numbers';
end

if nargin < 3
    labels = [];
end

nstates = size(states,1);
% Builds the binary matrix of the inclusion relation...
inclusion_relation = zeros(nstates);
for i = 1:nstates
    for j = 1:nstates
        inclusion_relation(i,j) = isequal(min(states(i,:),states(j,:)),states(i,:));   
    end
end
... and then the adjacency matrix of the covering relation
adjacency_matrix = inclusion_relation;
for i = 1:nstates
    for j = 1:nstates
        if inclusion_relation(i,j) == true
            for k = 1:nstates
                if k ~= i && k ~= j && inclusion_relation(i,k) == true &&...
                        inclusion_relation(k,j) == true
                    adjacency_matrix(i,j) = false;
                end
            end
        end
    end
end
% Construction of the Hasse diagram
G = digraph(adjacency_matrix,'OmitSelfLoops');
%card=sum(states,2);
hGraphPlot = plot(G,'k');
layout(hGraphPlot,'layered');
hGraphPlot.MarkerSize = 4;
hGraphPlot.ShowArrows = 'off';
set(gca,'XTick',[])
set(gca,'YTick',[])

if size(states,1) <= 100
    switch labelstyle
        case 'vectors'
            for i = 1:size(states,1)
                s = num2str(states(i,:));
                s(isspace(s))=[];
                labelnode(hGraphPlot,i,s);
            end
            
        case 'literals'            
            for i = 1:size(states,1)
                K = find(states(i,:));
                s = '{';
                for j = 1:length(K)
                    s = [s,char(96+K(j))];
                    if j < length(K)
                        s = [s,','];
                    end
                end
                s = [s,'}'];
                labelnode(hGraphPlot,i,s);
            end
        case 'numbers'
            for i = 1:size(states,1)
                K = find(states(i,:));
                s = '{';
                for j = 1:length(K)
                    s = [s,num2str(K(j))];
                    if j < length(K)
                        s = [s,','];
                    end
                end
                s = [s,'}'];
                labelnode(hGraphPlot,i,s);
            end 
        case 'custom'
            for i = 1:size(states,1)
                K = find(states(i,:));
                s = '{';
                for j = 1:length(K)
                    s = [s,labels{K(j)}];
                    if j < length(K)
                        s = [s,','];
                    end
                end
                s = [s,'}'];
                labelnode(hGraphPlot,i,s);
            end

        case 'cardinality'
            for i = 1:size(states,1)
                card=sum(states(i,:));
                labelnode(hGraphPlot,i,card);
            end
    end
end