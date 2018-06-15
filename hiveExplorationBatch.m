% Test Script for exploring Hives, uses parallel sampling
% dim must be >2

close all
clear
clc

%% Vector space dimension, can be an array (eg 5:5:20)
m = 4;

%% Total number of samples for each dimension
sampleMax = 1000;

%% Tollerance
tol = 10^-8;

%% Establish storage
hiveCheck = zeros(1,length(m));

%% Flags must be unset%
flagData.waitbar = 0; %
flagData.verb = 0;    %
flagData.parallel = 0;%
%%%%%%%%%%%%%%%%%%%%%%%

%% Setup waitbar and counters
h = waitbar(0,'Progress'); 
counter = 1;
    
%% Loop over dimensions
for n = m
    
    hiveProb = 0;

    %% Parallel loop over samples
    tests = 0;
    for q = 1:sampleMax

        
        %% Setup matrix ensembles
        
        %{
        L = GOEGenerator(n);
        M = GOEGenerator(n);
        %}
        
        %{
        %Normal SPD
        L = randn(n);
        L = L'*L;
        M = randn(n);
        M = M'*M;
        %}   
           
        %{
        M = diag(sort(randi(50,[1,n]),'descend'));
        L = sort(randi(50,[1,n]),'descend');
        sw = L(1);
        L(1) = L(end);
        L(end) = sw;
        L = diag(L);
        
        %}
        
        %Commuting Hermitian Matrices
        %{
        O = gallery('orthog',n,2);
        L = O' * diag(rand(n,1)) * O;
        M = O' * diag(rand(n,1)) * O;
        %}
        
        %Hermitian Matrices positive spectra
        %{
        O = gallery('orthog',n,2);
        L = O' * diag(randi(50,[n,1])) * O;
        O = gallery('orthog',n,1);
        M = O' * diag(randi(50,[n,1])) * O;
        %}
        
        %% Perform optimization with up to 5 re-runs to rule out local mins or instabilities
        
        %%{
        hiveLogical = false;
        
        trials = 0;
        while ~hiveLogical && trials < 100

            Hijk = AWHiveParallel(L,M,flagData);
            
            if ~isempty(Hijk)
                
                hiveLogical = rhombusCheckF(Hijk,n);
                
            end
            %%{
            trials = trials + 1;
            
            if trials == 100
                
                if isempty(Hijk)
                    
                    warning('empty')
                    
                elseif ~hiveLogical
                    
                    warning('failure')
                    
                end
                
            end
            %}
            tests = tests + 1;
        end
        

        hiveProb = hiveProb + hiveLogical;
        
        %}

    end    
    
    hiveCheck(counter) = hiveProb;
    counter = counter + 1;
    waitbar(counter/length(m))
    
end

%% Close waitbar
close(h)

%% Compute means and confidence intervals
hiveCheck = hiveCheck./sampleMax;
hiveCI = sqrt(hiveCheck.*(1-hiveCheck))./(sqrt(sampleMax)) * 1.96;

%% Plots

%Plot probability of forming a hive
figure(1)
errorbar(m,hiveCheck,hiveCI,'o-')
title(sprintf('Hive Probability for Independent Commuting Real Symmetric Matrices for %d Samples \n and 95%% Confidence Intervals',sampleMax))
xlabel('Matrix Dimension')
ylabel('Probability')