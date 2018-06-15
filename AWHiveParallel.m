function Hijk = AWHiveParallel(L,M,flagData)
%% Computes the Hijk hive coefficients
%based on the Appleby-Whitehead construction given an input of a pair of
%nxn real symmetric matrices L and M.
%Hijk is reported as a sparse matrix with rows corresponding to k+1, collumns i+1,
%such that L is restricted to a subspace of dimension i+k while M is
%restricted to subspace k; j = n-i-k is implicit.

%% Set some shortcuts and storage arrays for sparse indexing and sparse values
n = length(L);

%% Compute and store Eigenvalues of L, M, and L+M in Parallel for speed if possible

if flagData.parallel
    
    eVals = zeros(3,n);
    
    parfor i = 1:3

        switch i

            case 1

                eVals(i,:) = flip(eig(L))';

            case 2

                eVals(i,:) = flip(eig(M))';
                
            case 3
                
                eVals(i,:) = flip(eig(L+M))';

        end

    end 
    
    eVals = eVals';

else
    
    eVals = zeros(n,3);
    eVals(:,1) = flip(eig(L));
    eVals(:,2) = flip(eig(M));
    eVals(:,3) = flip(eig(L+M));
    
end

eVals = sort(eVals,'descend');
%% k = 0, i > 0: 
% The sum of the largest i eigenvalues gives the Hij0:

HijkV = [zeros(n,1), (1:n)', cumsum(eVals(:,1))];

%% k = n-i: L is always full rank,
% Get sum of evals of M from smallest space to largest (i goes from 1 to n) plus trace of L; 

HijkV = [HijkV; [(n-(0:(n-1)))',(0:(n-1))',sum(eVals(:,1)) + flip(cumsum(eVals(:,2)))]];

%% i=0, 0<k<n: 
% Mutually projected into same subspace, sum of evals of L+M with varying k

temp = cumsum(eVals(:,3));
HijkV= [HijkV;[(1:(n-1))',zeros(n-1,1), temp(1:(end-1))]];

%% Full product-grassmann trust-region for the rest, in parallel if possible

if flagData.parallel
   
    verb = flagData.verb;
    poolObj = gcp;
    
    for k = 1:n
        
        f(k) = parfeval(poolObj,@parHiveFctn,1,L,M,n,k,verb);
        
    end
    
    % This structure is used to abort the parallel calculations in the
    % event that one constant has an optimization that doesn't converge
    for k = 1:n
       
        [~,value] = fetchNext(f);
        
        if any(isinf(value(:,3)))
            
            cancel(f)
            warning('Optimization was non-convergent')
            Hijk = [];
            return
            
        else
            
            HijkV = [HijkV;value];
            
        end
        
    end
    
else
    
    %Waitbar initialize
    if flagData.waitbar

        h = waitbar(0,'Progress');

    end

    currentIndex = 3*n+1;
    
    for k = 1:n

        for i = 1:(n-k)

            %Skip, already done above
            if k == n-i

                continue

            end

            tmp = findHive(L,M,k,i,flagData.verb);
            
            if tmp == inf
                
                warning('Optimization was non-convergent')
                Hijk = [];
                return
                
            else
            
                HijkV(currentIndex,:) = [k,i,tmp];
                
            end
            
            currentIndex = currentIndex + 1;

        end
        
        if flagData.waitbar

            waitbar(k/n);

        end

    end
    
    if flagData.waitbar

        close(h);

    end
    
end

%% Assign into sparse matrix
Hijk = sparse(HijkV(:,1)+1,HijkV(:,2)+1,HijkV(:,3));

end

function out = parHiveFctn(L,M,n,k,verb)

out = zeros((n-k),3);

for i = 1:(n-k)

    %Skip, already done above
    if k == n-i

        continue

    end

    out(i,:) = [k,i,findHive(L,M,k,i,verb)];

end  

end
