function [pat,freq,niter]=patfreq(sample)
pat = [];
freq = [];
niter = 0;
visited = false(size(sample,1),1);
for i = 1:size(sample,1)    
    if ~visited(i)     
        f = 0;
        for j = i:size(sample,1)  
            if isequal(sample(i,:),sample(j,:))
                f = f+1;
                visited(j) = true;                
            end            
        end    
        pat = [pat;sample(i,:)];
        freq = [freq;f];
        
    end   
end