function out = findHive0(L,M,k,i)
%% Helper function for computing the proposed hive coefficient Hijk with either i or k = 0
% Called by the findHive general function.

global options

n = length(L);

%% Create the manifold structure.
if i == 0
    
    N = L+M;
    m = grassmannfactory(n,k);
    
else
    
    N = L;
    m = grassmannfactory(n,i);
    
end

%% Define the problem structure--manifold, cost function, and Euclidean gradient
problem.M = m;
problem.cost = @cost;
problem.egrad = @egrad;
problem.ehess = @ehess;

%% Define the problem cost function and its Euclidean gradient.

%Cost function
function [f,store] = cost(x,store)

    if ~isfield(store,'X')
        store.X = inv(x'*x);
    end
    if ~isfield(store, 'proj')
        store.proj =x*store.X*x';
    end

    f  = -trace(store.proj*N*store.proj);

end

%Euclidean gradient
function [g,store] = egrad(x,store)

    if ~isfield(store,'X')
        store.X = inv(x'*x);
    end
    if ~isfield(store, 'proj')
        store.proj = x*store.X*x';
    end

    g = -2*(eye(n)-store.proj)*N*x*store.X;

end

%Euclidean Hessian
function [h,store] = ehess(x,xd,store)

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


%% Solve.
[~, out, info, options] = trustregions(problem,[],options);

redoCounter = 0;
    
while  numel([info.gradnorm]) > 50 &&  info(end).gradnorm==info(end-10).gradnorm && redoCounter < 10

    if options.verbosity > 0
        disp('Maxed out--redo')
    end
    [~, out, info, options] = trustregions(problem,[],options);
    redoCounter = redoCounter + 1;

end

if numel(info) > 50 && info(end).gradnorm==info(end-10).gradnorm
    
    out = inf;
    
else
    
    out = -out;
    
end

end