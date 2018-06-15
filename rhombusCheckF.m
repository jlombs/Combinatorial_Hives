function out = rhombusCheckF(Hijk,n)
%% Function for checking rhombus inequalities for AWHives,
%outputs logical if it forms a good hive (1-hive, 0 fail).
%Fast logical version of rhombusCheck for use with serial processing.
%n is the vector space dimension, and Hijk is a sparse hive matrix.

out = 1;

%% Loop over all possible indices and check the 3 rhombuses at each point

for k = 1:(n+1)
    
    for i = 1:(n+2-k)
        
        %% Right Inequality
        
            if i <= (n-k) && Hijk(k,i) + Hijk(k+1,i+1) - Hijk(k,i+1) - Hijk(k+1,i) > 0
                    
                out = 0;
                return
                
            end

        %% Top Inequality
        
            if  i > 2 && Hijk(k,i) + Hijk(k+1,i-2) - Hijk(k,i-1) - Hijk(k+1,i-1) > 0
                    
                out = 0;
                return
                
            end

        %% Left Inequality
        
            if k > 2 && Hijk(k,i) + Hijk(k-2,i+1) - Hijk(k-1,i) - Hijk(k-1,i+1) > 0
                    
                out = 0;
                return
                
            end
            
    end
  
end

end
