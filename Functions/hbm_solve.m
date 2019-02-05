function sol = hbm_solve(hbm,problem,w0,A,x0)
NDof = problem.NDof;
Nfft  = hbm.harm.Nfft;

if nargin < 5 || isempty(x0)
    x0 = zeros(hbm.harm.NComp*NDof);
end
x0 = packdof(x0,hbm.harm.iRetain);

%setup the problem for IPOPT
%we have an initial guess vector
U = A*feval(problem.excite,hbm,problem,w0*hbm.harm.rFreqRatio);

hbm.bIncludeNL = true;

u = packdof(U);

iRetain = hbm.harm.iRetain;
NComp = hbm.harm.NComp;
Nhbm  = hbm.harm.NRetain;

if isfield(problem,'xhscale')
    xdc  = problem.x0scale(:)';
    xmax = problem.xhscale(:)';
    xscale = [xdc; repmat(xmax,hbm.harm.NFreq-1,1)*(1+1i)];
    problem.xscale = packdof(xscale,hbm.harm.iRetain)*sqrt(length(xscale));
    problem.Fscale = problem.xscale*0+1;
    problem.Jscale = (1./problem.Fscale(:))*problem.xscale(:)';
else
    problem = hbm_scaling(problem,hbm,x0);
end

%and actually solve it
hbm.max_iter = 3;
bSuccess = false;

constr_tol = 1E-6;
maxit = 20;

X0 = x0./problem.xscale;
% F0 = feval(funcs.constraints,X0,options.auxdata);
% J0 = feval(funcs.jacobian,X0,options.auxdata);

attempts = 0;
while ~bSuccess && attempts < hbm.max_iter
    switch hbm.options.solver
        case 'fsolve'
             fun_constr = @(x)hbm_constraints(x,hbm,problem,w0,u);
             options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',true,'FunctionTolerance',constr_tol,'MaxIterations',maxit);
             [X,~,EXITFLAG,OUTPUT] = fsolve(fun_constr,X0,options);
             bSuccess = EXITFLAG == 1;
             iter = OUTPUT.iterations + 1;
        case 'ipopt'
            options.jacobian = @hbm_jacobian;
            options.jacobianstructure = hbm.sparsity;
            options.print_level = 5;
            options.max_iter = maxit;           
            options.constr_viol_tol = constr_tol;
            [X, info] = fipopt('',X0,@hbm_constraints,options,hbm,problem,w0,u);
            bSuccess = any(info.status == [0 1]);
            iter = info.iter;
    end
    X0 = X + 1E-8*rand(Nhbm,1);
    attempts = attempts + 1;
end
if ~bSuccess
    X = X + NaN;
end
% F = feval(funcs.constraints,X,options.auxdata);

x = X.*problem.xscale;
sol.X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof,iRetain);
sol.w0 = w0;

w = w0*hbm.harm.rFreqRatio;
sol.w = w;

f = hbm_output3d(hbm,problem,w,u,x);
sol.F = unpackdof(f,hbm.harm.NFreq-1,problem.NOutput);

sol.U = U;
sol.A = A;

%floquet multipliers
sol.L = hbm_floquet(hbm,problem,w,u,x);

sol.it = iter;

function [c,J] = hbm_constraints(X,hbm,problem,w0,u)
%unpack the inputs
x = X.*problem.xscale;
w = w0*hbm.harm.rFreqRatio;
c = hbm_balance3d('func',hbm,problem,w,u,x);
c = c ./ problem.Fscale;
if nargout > 1
    J = hbm_jacobian(X,hbm,problem,w0,u);
end

function J = hbm_jacobian(X,hbm,problem,w0,u)
%unpack the inputs
% J = jacobian(@hbm_constraints,x,auxdata);
% J = sparse(J);
% return;

x = X.*problem.xscale;
w0 = w0*hbm.harm.rFreqRatio;
J = hbm_balance3d('jacob',hbm,problem,w0,u,x);
J = J .* problem.Jscale;