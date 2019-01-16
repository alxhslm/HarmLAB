function sol = hbm_solve(hbm,problem,w0,A,x0)
NDof = problem.NDof;
NAlg = problem.NAlg;
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
iter = 0;
hbm.max_iter = 3;
bSuccess = false;

X0 = x0./problem.xscale;
% F0 = feval(funcs.constraints,X0,options.auxdata);
% J0 = feval(funcs.jacobian,X0,options.auxdata);

while ~bSuccess && iter < hbm.max_iter
    switch hbm.options.solver
        case 'nleqn'
             [X,info] = nleqn(@(x)hbm_constraints(x,hbm,problem,w0,u),X0);
             bSuccess = ~info;
        case 'fsolve'
             fun_constr = @(x)hbm_constraints(x,hbm,problem,w0,u);
             options = optimoptions('fsolve','Display','iter');
             [X,~,EXITFLAG] = fsolve(fun_constr,X0,options);
             bSuccess = EXITFLAG == 1;
        case 'ipopt'
            options.jacobian = @hbm_jacobian;
            options.jacobianstructure = hbm.sparsity;
            options.print_level = 5;
            options.max_iter = 15;           
            
            [X, info] = fipopt('',X0,@hbm_constraints,options,hbm,problem,w0,u);
            bSuccess = any(info.status == [0 1]);
    end
    X0 = X + 1E-8*rand(Nhbm+prod(Nfft)*NAlg,1);
    iter = iter + 1;
end
if ~bSuccess
    X = X + NaN;
end
% F = feval(funcs.constraints,X,options.auxdata);

x = X.*problem.xscale;
sol.X = unpackdof(x(1:Nhbm),hbm.harm.NFreq-1,problem.NDof,iRetain);
sol.xalg = permute(reshape(x(Nhbm+1:end,:),NAlg,Nfft(1),Nfft(2)),[2 3 1]);
sol.w0 = w0;

w = w0*hbm.harm.rFreqRatio;
sol.w = w;

f = hbm_output3d(hbm,problem,w,u,x);
sol.F = unpackdof(f,hbm.harm.NFreq-1,problem.NOutput);

sol.U = U;
sol.A = A;

%floquet multipliers
sol.L = floquetMultipliers(hbm,problem,w,u,x);

function c = hbm_constraints(X,hbm,problem,w0,u)
%unpack the inputs
x = X.*problem.xscale;
w0 = w0*hbm.harm.rFreqRatio;
c = hbm_balance3d('func',hbm,problem,w0,u,x);
c = c ./ problem.Fscale;

function J = hbm_jacobian(X,hbm,problem,w0,u)
%unpack the inputs
% J = jacobian(@hbm_constraints,x,auxdata);
% J = sparse(J);
% return;

x = X.*problem.xscale;
w0 = w0*hbm.harm.rFreqRatio;
J = hbm_balance3d('jacob',hbm,problem,w0,u,x);
J = J .* problem.Jscale;