%% Test Script for checking agreement of boundaries with optimization

close all
clc
clear

%% Set input matrices and optmization flags
m = 8;

flagData.waitbar = 0;
flagData.verb = 0;
flagData.parallel = 1;

L = GOEGenerator(m);
M = GOEGenerator(m);

if flagData.verb > 0
    
    fprintf('\n Optimization Feedback \n')
    
end

%% Compute percent abs error in each term

% Compute the known true values
v = zeros(3*m-1,1);
eVals = zeros(m,3);
eVals(:,1) = flip(eig(L));
eVals(:,2) = flip(eig(M));
eVals(:,3) = flip(eig(L+M));
eVals = sort(eVals,'descend');

v(1:m) = cumsum(eVals(:,1));
v((m+1):(2*m)) = sum(eVals(:,1)) + flip(cumsum(eVals(:,2)));
temp = cumsum(eVals(:,3));
v((2*m+1):end) = temp(1:(end-1));
    
%Compute the opt versions by reference index
optVals = zeros(3*m-1,1);
ki = [zeros(m,1), (1:m)';(m-(0:(m-1)))',(0:(m-1))';(1:(m-1))',zeros(m-1,1)];
parfor p = 1:(3*m-1)
    
    optVals(p) = findHive(L,M,ki(p,1),ki(p,2),1);
    
end

%Compute errors and CI
errors = 100.*(v - optVals)./v;

stdev = std(errors)*1.96;

%% Plot Results
figure(1)
scatter(1:length(errors),errors,'*r')
hold on
lineX = 0:.1:(length(errors)+1);
plot(lineX,stdev*(ones(1,length(lineX))),'-k')
plot(lineX,-stdev*(ones(1,length(lineX))),'-k')
legend([{'Errors'},{'95% UCI'},{'95% LCI'}])
xlabel('Element Number')
ylabel('Percent Error (%)')
title('Percent Error Between Analytical Boundary Values and Optimization Derived Values')