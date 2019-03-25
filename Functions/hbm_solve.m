function sol = hbm_solve(hbm,problem,w,A,X0)
problem.type = 'solve';

NDof = problem.NDof;

if nargin < 5 || isempty(X0)
    X0 = zeros(hbm.harm.NComp*NDof);
end
if ~isvector(X0)
    x0 = packdof(X0,hbm.harm.iRetain);
else
    x0 = X0;
end

%setup the problem for IPOPT
%we have an initial guess vector
w0 = w*hbm.harm.rFreqRatio;
U = A*feval(problem.excite,hbm,problem,w0);

hbm.bIncludeNL = true;

u = packdof(U);

iRetain = hbm.harm.iRetain;
NComp = hbm.harm.NComp;
Nhbm  = hbm.harm.NRetain;

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

z0 = x0;
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
    Z0 = Z + 1E-8*rand(Nhbm,1);
    attempts = attempts + 1;
end
if ~bSuccess
    Z = Z + NaN;
end
z = Z.*problem.Zscale;

x = z;

sol.w = w;
sol.A = A;
sol.X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof,iRetain);
sol.U = U;
sol.F = hbm_output3d(hbm,problem,w0,sol.U,sol.X);

%floquet multipliers
sol.L = hbm_floquet(hbm,problem,w0,u,x);

sol.it = iter;

function [c,J] = hbm_constraints(Z,hbm,problem,w0,u)
%unpack the inputs
x = Z.*problem.xscale;
w = w0*hbm.harm.rFreqRatio;
c = hbm_balance3d('func',hbm,problem,w,u,x);
c = c ./ problem.Fscale;
if nargout > 1
    J = hbm_jacobian(Z,hbm,problem,w0,u);
end

function J = hbm_jacobian(Z,hbm,problem,w0,u)
x = Z.*problem.xscale;
w0 = w0*hbm.harm.rFreqRatio;
J = hbm_balance3d('jacob',hbm,problem,w0,u,x);
J = J .* problem.Jscale;