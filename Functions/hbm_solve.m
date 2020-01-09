function sol = hbm_solve(hbm,problem,w,A,X0)
problem.type = 'solve';

if nargin < 5 || isempty(X0)
    X0 = zeros(hbm.harm.NFreq,problem.NDof);
end
if isvector(X0) 
    if length(X0) == hbm.harm.NComp*problem.NDof
        x0 = X0;
        X0 = unpackdof(x0,hbm.harm.NHarm,problem.NDof);
    else
        error('Wrong size for X0');
    end
else
    if size(X0,1) == hbm.harm.NFreq && size(X0,2) == problem.NDof
        x0 = packdof(X0);
    else
        error('Wrong size for X0');
    end
end

%setup the problem for IPOPT
%we have an initial guess vector
w0 = w * hbm.harm.rFreqRatio + hbm.harm.wFreq0;
U = A*feval(problem.excite,hbm,problem,w0);

hbm.bIncludeNL = true;

u = packdof(U);

NRetainNL  = hbm.harm.NRetainNL;

init.X = X0;
init.w = w;
init.A = A;

if isfield(problem,'xhscale')
    xdc  = problem.x0scale(:)';
    xmax = problem.xhscale(:)';
    xscale = [xdc; repmat(xmax,hbm.harm.NFreq-1,1)*(1+1i)];
    problem.xscale = packdof(xscale,hbm.harm.iRetain)*sqrt(length(xscale));
    problem.Fscale = problem.xscale*0+1;
    problem.Jscale = (1./problem.Fscale(:))*problem.xscale(:)';
else
    problem = hbm_scaling(problem,hbm,init);
end

%and actually solve it
hbm.max_iter = 3;
bSuccess = false;

constr_tol = 1E-6;
maxit = 20;

z0 = x0(hbm.harm.iNL);
Z0 = z0./problem.xscale;

attempts = 0;
while ~bSuccess && attempts < hbm.max_iter
    switch hbm.options.solver
        case 'fsolve'
             fun_constr = @(x)hbm_constraints(x,hbm,problem,w,u);
             options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',true,'FunctionTolerance',constr_tol,'MaxIterations',maxit);
             [Z,~,EXITFLAG,OUTPUT] = fsolve(fun_constr,Z0,options);
             bSuccess = EXITFLAG == 1;
             iter = OUTPUT.iterations + 1;
        case 'ipopt'
            options.jacobian = @hbm_jacobian;
            options.jacobianstructure = hbm.sparsity;
            options.print_level = 5;
            options.max_iter = maxit;           
            options.constr_viol_tol = constr_tol;
            [Z, info] = fipopt('',Z0,@hbm_constraints,options,hbm,problem,w,u);
            bSuccess = any(info.status == [0 1]);
            iter = info.iter;
    end
    Z0 = Z + 1E-8*rand(NRetainNL,1);
    attempts = attempts + 1;
end
if ~bSuccess
    Z = Z + NaN;
end
z = Z.*problem.Zscale;

x = hbm_recover3d(hbm,problem,w,u,z);
X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof,hbm.harm.iRetain);

sol.w = w;
sol.A = A;
sol.X = X;
sol.U = U;
sol.F = hbm_output3d(hbm,problem,w,sol.U,sol.X);

% floquet multipliers
sol.L = hbm_floquet(hbm,problem,w,sol.U,sol.X);

sol.it = iter;

function [c,J] = hbm_constraints(Z,hbm,problem,w,u)
%unpack the inputs
x = Z.*problem.xscale;
c = hbm_balance3d('func',hbm,problem,w,u,x);
c = c ./ problem.Fscale;
if nargout > 1
    J = hbm_jacobian(Z,hbm,problem,w,u);
end

function J = hbm_jacobian(Z,hbm,problem,w,u)
x = Z.*problem.xscale;
J = hbm_balance3d('jacob',hbm,problem,w,u,x);
J = J .* problem.Jscale;