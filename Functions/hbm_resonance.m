function sol = hbm_resonance(hbm,problem,w0,A,xIn)
NDof = problem.NDof;
NAlg = problem.NAlg;
Nfft  = hbm.harm.Nfft;

%setup the problem for IPOPT

hbm.bIncludeNL = 1;
NComp = hbm.harm.NComp;
Nhbm  = hbm.harm.NRetain;

% %first solve @ w0
% x0 = hbm_solve(hbm,problem,w0,A,xIn);
x0 = xIn;
x0 = packdof(x0);

if isfield(problem,'xhscale')
    xdc  = problem.x0scale(:)';
    xmax = problem.xhscale(:)';
    xscale = [xdc; repmat(xmax,hbm.harm.NFreq-1,1)*(1+1i)];
    problem.xscale = packdof(xscale,hbm.harm.iRetain);
    problem.Fscale = problem.xscale*0+1;
    problem.wscale = w0;
    problem.Xscale = [problem.xscale; problem.wscale];
else
    problem = hbm_scaling(problem,hbm,x0,w0);
end

problem.Jscale = (1./problem.Fscale(:))*problem.Xscale(:)';

X0 = [x0; w0]./problem.Xscale;

%and actually solve it
iter = 0;
hbm.max_iter = 4;
bSuccess = false;
Jstr = [hbm.sparsity ones(NDof*NComp,1)];

% obj0 = hbm_obj(X0,hbm,problem,A);
% G0 = hbm_grad(X0,hbm,problem,A);
% J0 = hbm_jacobian(X0,hbm,problem,A);
% F0 = hbm_constr(X0,hbm,problem,A);
% 
% dHdw0 = G0(end) - G0(1:end-1)*(J0(:,1:end-1)\J0(:,end));

constr_tol = 1E-6;
opt_tol = 1E-6;
maxit = 20;

while ~bSuccess && iter < hbm.max_iter
    switch hbm.options.solver
        case 'fsolve'
            fun_obj = @(x)hbm_fsolve_obj(x,hbm,problem,A);
            fun_constr = @(x)hbm_fsolve_constr(x,hbm,problem,A);
            options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'Display','iter',...
                'OptimalityTolerance',opt_tol,'ConstraintTolerance',constr_tol,'MaxIterations',maxit);
            [X,~,EXITFLAG] = fmincon(fun_obj,X0,[],[],[],[],[],[],fun_constr,options);
            bSuccess = EXITFLAG == 1;
        case 'ipopt'
            options.jacobianstructure  = Jstr;
            options.jacobian = @hbm_jacobian;
            options.gradient = @hbm_grad;
            options.print_level = 5;
            options.max_iter = maxit;
            options.tol = opt_tol;
            options.constr_viol_tol = constr_tol;
            [X, info] = fipopt(@hbm_obj,X0,@hbm_constr,options,hbm,problem,A);           
            bSuccess = any(info.status == [0 1]);
    end
    X0 = X+rand(length(X),1)*1E-8;
    iter = iter + 1;
end
if ~bSuccess
    X = X + NaN;
end

% obj = hbm_obj(X,hbm,problem,A);
% G = hbm_grad(X,hbm,problem,A);
% J = hbm_jacobian(X,hbm,problem,A);
% F = hbm_constr(X,hbm,problem,A);
% 
% dHdw = G(end) - G(1:end-1)*(J(:,1:end-1)\J(:,end));

w0 = abs(X(end).*problem.wscale);
w = w0*hbm.harm.rFreqRatio;
x = X(1:end-1).*problem.xscale;

sol.X = unpackdof(x(1:Nhbm),hbm.harm.NFreq-1,NDof,hbm.harm.iRetain);
sol.xalg = squeeze(permute(reshape(x(Nhbm+1:end,:),NAlg,Nfft(1),Nfft(2)),[2 3 1]));

sol.w0 = w0;
sol.w = w;

sol.U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(sol.U);

sol.A = A;

f = hbm_output3d(hbm,problem,w,u,x);
sol.F = unpackdof(f,hbm.harm.NFreq-1,problem.NOutput);

%floquet multipliers
sol.L = floquetMultipliers(hbm,problem,w,u,x);

function obj = hbm_obj(X,hbm,problem,A)
x = X(1:end-1).*problem.xscale;
w = X(end).*problem.wscale;
w0 = w * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

H = hbm_objective('func',hbm,problem,w0,x,u);
obj = - problem.res.sign * H;

function G = hbm_grad(X,hbm,problem,A)
x = X(1:end-1).*problem.xscale;
w = X(end).*problem.wscale;
w0 = w * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

[Dx, Dw] = hbm_objective({'jacobX','derivW'},hbm,problem,w0,x,u);
G = -problem.res.sign*[Dx Dw];

G = G.*problem.Xscale(:)';

function c = hbm_constr(X,hbm,problem,A)
x = X(1:end-1).*problem.xscale;
w = X(end).*problem.wscale;
w0 = w * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

c = hbm_balance3d('func',hbm,problem,w0,u,x);

function J = hbm_jacobian(X,hbm,problem,A)
x = X(1:end-1).*problem.xscale;
w = X(end).*problem.wscale;
w0 = w * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

Jx_nl = hbm_balance3d('jacob',hbm,problem,w0,u,x);
Dw_nl = hbm_balance3d('derivW',hbm,problem,w0,u,x);


J = [Jx_nl  Dw_nl];
J = J .* problem.Jscale;

function [obj,G] = hbm_fsolve_obj(X,hbm,problem,A)
obj = hbm_obj(X,hbm,problem,A);
if nargout > 1
    G =  hbm_grad(X,hbm,problem,A)';
end

function [c,ceq,Jc,Jceq] = hbm_fsolve_constr(X,hbm,problem,A)
ceq = hbm_constr(X,hbm,problem,A);
c = [];
% c = ceq(end);
% ceq = ceq(1:end-1);
if nargout > 2
    Jceq = hbm_jacobian(X,hbm,problem,A);
    Jc = zeros(0,length(X));
    
    Jc = Jc';
    Jceq = Jceq';
%     Jc = Jceq(end,:)';
%     Jceq = Jceq(1:end-1,:)';
end

