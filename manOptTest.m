function out = manOptTest()
% Script to test the accuracy of the analytical gradient and hessian from the Grassmann map
% A variety of plots can be uncommented at the end, for examining the
% optimization in detail

close all
clear

%% Setup matrix dimensions and ensembles

n = 15;

L = GOEGenerator(n);
M = L;


%% Select subspace dimensions to test
% Based on our notation, for max(U,V) in F_k,(i+k) of tr(L|V) + Tr(M|U)

k = 8;
i = 6;


%% Check inputs
if i + k > n || i < 0 || k < 0 || k == n || i == n
    
    error('Invalid hive indices')
    
end

%% H0n0 is normalized to 0

if k == 0 && i == 0
    
    out = 0;
    return
    
end

%% Handle special cases of i == 0 or k == 0 where the manifold isn't a product

if i == 0 || k == 0

    %% Create the manifold structure.
    if i == 0

        N = L+M;
        m = grassmannfactory(n,k);
        manifold = m;

    else

        N = L;
        m = grassmannfactory(n,i);
        manifold = m;

    end

    %% Define the problem structure--manifold, cost function, and Euclidean gradient
    problem.M = m;
    problem.cost = @cost0;
    problem.egrad = @egrad0;
    problem.ehess = @ehess0;
    
else
    
    %% Create the manifold structure.
    tuple.mB = grassmannfactory(n,k);
    tuple.mAT = grassmannfactory(n,i);
    manifold = productmanifold(tuple);

    %% Define the problem structure--manifold, cost function, and Euclidean gradient
    problem.M = manifold;
    problem.cost = @cost;
    problem.egrad = @egrad;
    problem.ehess = @ehess;

end
    

%% Define the problem cost function and its Euclidean gradient for the special case

%Cost function
function [f,store] = cost0(x,store)

    if ~isfield(store,'X')
        store.X = inv(x'*x);
    end
    if ~isfield(store, 'proj')
        store.proj =x*store.X*x';
    end
    
    f  = -trace(store.proj*N*store.proj);

end

%Euclidean gradient
function [g,store] = egrad0(x,store)

    if ~isfield(store,'X')
        store.X = inv(x'*x);
    end
    if ~isfield(store, 'proj')
        store.proj = x*store.X*x';
    end

    g = -2*(eye(n)-store.proj)*N*x*store.X;

end

%Euclidean Hessian
function [h,store] = ehess0(x,xd,store)

    if ~isfield(store,'X')
        store.X = inv(x'*x);
    end
    if ~isfield(store, 'proj')
        store.proj =x*store.X*x';
    end
 
    Xd = -store.X*(xd'*x+x'*xd)*store.X;
    projd = xd*store.X*x'+x*(Xd*x'+ store.X*xd');

    h = -2*(-projd*N*x*store.X+(eye(n)-store.proj)*N*(xd*store.X+x*Xd));

end

%% Define the problem cost function and its Euclidean gradient for the i,k > 0 case

%Cost function
function [f,store] = cost(x,store)

    if ~isfield(store,'A')
        store.A = [x.mB,x.mAT];
    end
    if ~isfield(store,'B')
        store.B = x.mB;
    end
    if ~isfield(store,'XA')
        store.XA = inv(store.A'*store.A);
    end
    if ~isfield(store,'XB')
        store.XB = inv(store.B'*store.B);
    end
    if ~isfield(store, 'projV')
        store.projV = store.A*store.XA*store.A';
    end
    if ~isfield(store, 'projU')
        store.projU = store.B*store.XB*store.B';
    end

    f  = -(trace(store.projV*L*store.projV)+trace(store.projU*M*store.projU));

end

%Euclidean gradient
function [g,store] = egrad(x,store)

    if ~isfield(store,'A')
        store.A = [x.mB,x.mAT];
    end
    if ~isfield(store,'B')
        store.B = x.mB;
    end
    if ~isfield(store,'XA')
        store.XA = inv(store.A'*store.A);
    end
    if ~isfield(store,'XB')
        store.XB = inv(store.B'*store.B);
    end
    if ~isfield(store, 'projV')
        store.projV = store.A*store.XA*store.A';
    end
    if ~isfield(store, 'projU')
        store.projU = store.B*store.XB*store.B';
    end
    if ~isfield(store,'X13')
        store.X13 = padarray(eye(k),[i,0],'post');
    end
    if ~isfield(store,'X24')
        store.X24 = padarray(eye(i),[k,0],'pre');
    end
    if ~isfield(store,'LTilde')
        store.LTilde = L*store.projV;
    end
    if ~isfield(store,'gradTemp')
        store.gradTemp = (eye(n)-store.projV)*(store.LTilde+store.LTilde')*store.A;
    end

    g.mB = -2*((eye(n)-store.projU)*M*store.B*store.XB + store.gradTemp*store.XA*store.X13);
    g.mAT = -2*store.gradTemp*store.XA*store.X24;

end

%Euclidean Hessian
function [h,store] = ehess(x,xd,store)

    if ~isfield(store,'A')
        store.A = [x.mB,x.mAT];
    end
    if ~isfield(store,'B')
        store.B = x.mB;
    end
    if ~isfield(store,'XA')
        store.XA = inv(store.A'*store.A);
    end
    if ~isfield(store,'XB')
        store.XB = inv(store.B'*store.B);
    end
    if ~isfield(store, 'projV')
        store.projV = store.A*store.XA*store.A';
    end
    if ~isfield(store, 'projU')
        store.projU = store.B*store.XB*store.B';
    end
    if ~isfield(store,'X13')
        store.X13 = padarray(eye(k),[i,0],'post');
    end
    if ~isfield(store,'X24')
        store.X24 = padarray(eye(i),[k,0],'pre');
    end
    if ~isfield(store,'LTilde')
        store.LTilde = L*store.projV;
    end
    if ~isfield(store,'gradTemp')
        store.gradTemp = (eye(n)-store.projV)*(store.LTilde+store.LTilde')*store.A;
    end
    
    Ad = [xd.mB,xd.mAT];
    Xd = -store.XA*(Ad'*store.A+store.A'*Ad)*store.XA;
    projVd = Ad*store.XA*store.A'+store.A*Xd*store.A'+store.A*store.XA*Ad';
    BInvd = -store.XB*(xd.mB'*store.B+store.B'*xd.mB)*store.XB;
    projUd = xd.mB*store.XB*store.B'+store.B*BInvd*store.B'+store.B*store.XB*xd.mB';

    LTilded = L*projVd;
    
    temp = (-projVd*(store.LTilde+store.LTilde')*store.A + (eye(n)-store.projV)*((LTilded+LTilded')*store.A + (store.LTilde+store.LTilde')*Ad))*store.XA...
        + store.gradTemp*Xd;

    h.mB = -2*(-projUd*M*store.B*store.XB+(eye(n)-store.projU)*M*(xd.mB*store.XB+store.B*BInvd)...
        + temp*store.X13);
    h.mAT = -2*temp*store.X24;
    
    %% Print this variable for to check for divergences
    %store.XA

end

%% Numerically check gradient consistency (optional--passed).
%{
figure()
checkgradient(problem)
figure()
checkhessian(problem)
%}

%% Solve in 2 ways for comparison
%%{
options.verbosity = 0;
options.stopfun = @mystopfun;
options.debug = 0;
options.userand = true;
[x, xcost, info, options] = trustregions(problem,[],options);

redoCounter = 0;

if numel(info) > 50
    
    while  info(end).gradnorm==info(end-10).gradnorm && redoCounter < 10

        if options.verbosity > 0
            disp('Maxed out--redo')
        end
        %{
        figure;
        semilogy([info.iter], [info.gradnorm], '.-');
        xlabel('Iteration number');
        ylabel('Norm of the gradient of f');
        %}
        [~, xcost, info, options] = trustregions(problem,[],options);
        redoCounter = redoCounter + 1;

    end

end

if numel(info) > 50 && info(end).gradnorm==info(end-10).gradnorm
    
    out = inf;
    
else
    
    out = -xcost;
    
end


% Display some statistics.
%{
figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');
%}

%% Plot pictures of cost function
%{
figure()
plotprofile(problem,manifold.rand(),manifold.rand(),-5:.05:5);
title('Cost Function along a Geodesic of the Manifold')
figure()
surfprofile(problem,x,manifold.rand(),manifold.rand(),-10:.3:10,-10:.3:10);
title('Cost Function as a Surface within a Neighborhood on the Manifold')
%}

%{
options.maxcostevals = 10000;

[~, xcost2, ~, ~] = pso(problem,[],options);



fprintf('Cost Function Optimum Difference between Trust Region and PSO Methods = %3.2e',abs(xcost-xcost2))
xcost
xcost2
%}
end

function stopnow = mystopfun(problem, x, info, last)
%% Optimization cut function
% If at least 50 iterations have elapsed and the current gradient norm is
% the same as it was 10 iterations ago, then this is characteristic of a
% divergence in our problem, and we can abort the optimization (or else the
% optimization would continue 'stuck' until timeout, which wastes
% computational time.

stopnow = (last >= 50 && info(last-10).gradnorm == info(last).gradnorm);
end