% Test Script for exploring Hives, single sample data focused
% dim m must be >2


close all
%clc
%clear

flagData.waitbar = 1;
flagData.verb = 0;
flagData.parallel = 1;

%{
% GOE matrices, sames
m = 10;
L = GOEGenerator(m);
M = L;
%}

%Hermitian Matrices
%{
m = 5;
O = gallery('orthog',m,2);
L = O' * diag(randi(10,[1,m])) * O;
M = O' * diag(randi(10,[1,m])) * O;
%}
        
%{
%GOE matrices, different
m = 6;
L = GOEGenerator(m);
M = GOEGenerator(m);
%}

%{
% Diagonal matrices with weakly decreasing positive integer values 
m = 4;
M = diag(sort(randi(10,[1,m]),'descend'));
L = diag(sort(randi(10,[1,m]),'ascend'));
%}

%{
% Diagonal matrices with weakly decreasing real values less than abs(x)
x = 500;
m = 20;
M = diag(sort(x*(2*rand([m,1])-1),'descend'));
L = diag(sort(x*(2*rand([m,1])-1),'descend'));
%}

%{
% Diagonal matrices with positive integer values less than x
m = 10;
x = 10;
M = diag(randi(x,[1,m]));
L = diag(randi(x,[1,m]));
%}

%{
% Diagonally Dominant SPD matrices
m = 10;
L = RSMGenerator(m);
L = L + m*eye(m);
M = RSMGenerator(m);
M = M + m*eye(m);
%}

%{
% Diagonally Dominant SPD matrices. sames
m = 20;
L = RSMGenerator(m);
L = L + m*eye(m);
M = L;
%}

%{
% Normally distributed SPD, sames
m = 20;
L = randn(m);
L = L'*L;
M = L;

%}

%{
% Normally distributed SPD
m = 20;
L = randn(m);
L = L'*L;
M = randn(m);
M = M'*M;

%}

%{
% Symmetric real matrices with integer eigenvalues <= x such that
% |mu|+|nu|=|lambda| but lambda are not necessarily integers...
m = 20;
x = 10;
V = orth(randn(m));
L = V*diag(randi(x,[m,1]))*V';
V = orth(randn(m));
M = V*diag(randi(x,[m,1]))*V';

%}

tol = 10^-6;

if flagData.verb > 0
    
    fprintf('\n Rough Optimization Feedback \n')
    
end

Hijk = AWHiveParallel(L,M,flagData);

%Hijk = round(Hijk);

%% Plot Hive and Rhombus Failures if Applicable
if ~isempty(Hijk)
    
    plotHive(Hijk,m,tol)
    
end