function out = findHive(L,M,k,i,verb)
%% Helper function for computing the proposed hive coefficient Hijk (for provided i,k)
%for real symmetric matrices L,M of dimension nxn based on the Appleby-Whitehead
%construction through continuum optimization on a product of grassmann
%manifolds. Uses manopt package.
%verb is the verbosity requested in the opt output

global options
    
options.verbosity = verb;
options.useRand = false;
options.stopfun = @mystopfun;

%% Normalize H0n0 to 0
if k == 0 && i == 0
    
    out = 0;
    return
    
end

%% Handle special cases of i == 0 or k == 0

if i == 0 || k == 0
    
    out = findHive0(L,M,k,i);
    return
    
end

%% k,i > 0

n = length(L);

%% Create the manifold structure.
tuple.mB = grassmannfactory(n,k);
tuple.mAT = grassmannfactory(n,i);
manifold = productmanifold(tuple);

%% Define the problem structure--manifold, cost function, and Euclidean gradient
problem.M = manifold;
problem.cost = @cost;
problem.egrad = @egrad;
problem.ehess = @ehess;

%% Define the problem cost function and its Euclidean gradient.

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

    f = -(trace(store.projV*L*store.projV)+trace(store.projU*M*store.projU));

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
    if ~isfield(store.shared,'X13')
        store.shared.X13 = padarray(eye(k),[i,0],'post');
    end
    if ~isfield(store.shared,'X24')
        store.shared.X24 = padarray(eye(i),[k,0],'pre');
    end
    if ~isfield(store,'LTilde')
        store.LTilde = L*store.projV;
    end
    if ~isfield(store,'gradTemp')
        store.gradTemp = (eye(n)-store.projV)*(store.LTilde+store.LTilde')*store.A;
    end

    g.mB = -2*((eye(n)-store.projU)*M*store.B*store.XB + store.gradTemp*store.XA*store.shared.X13);
    g.mAT = -2*store.gradTemp*store.XA*store.shared.X24;

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
    if ~isfield(store.shared,'X13')
        store.shared.X13 = padarray(eye(k),[i,0],'post');
    end
    if ~isfield(store.shared,'X24')
        store.shared.X24 = padarray(eye(i),[k,0],'pre');
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
        + temp*store.shared.X13);
    h.mAT = -2*temp*store.shared.X24;

end

%% Solve.

[~, out, info, options] = trustregions(problem,[],options);
%{
[~, fbest, ~, options] = pso(problem) ;
if fbest < out
    out = fbest;
end
%}
redoCounter = 0;
    
while  numel([info.gradnorm]) > 50 && info(end).gradnorm==info(end-10).gradnorm && redoCounter < 10

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

function stopnow = mystopfun(problem, x, info, last)
%% Optimization cut function
% If at least 50 iterations have elapsed and the current gradient norm is
% the same as it was 10 iterations ago, then this is characteristic of a
% divergence in our problem, and we can abort the optimization (or else the
% optimization would continue 'stuck' until timeout, which wastes
% computational time.

stopnow = (last >= 50 && info(last-10).gradnorm == info(last).gradnorm);
end