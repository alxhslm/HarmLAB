function sol = hbm_resonance(hbm,problem,w,A,X0)
problem.type = 'resonance';

NDof = problem.NDof;
Nfft  = hbm.harm.Nfft;

%setup the problem for IPOPT

hbm.bIncludeNL = 1;
NComp = hbm.harm.NComp;
Nhbm  = hbm.harm.NRetain;

%first solve @ w0
sol = hbm_solve(hbm,problem,w,A,X0);
X0 = sol.X;
x0 = packdof(X0);

init.X = X0;
init.w = w;
init.A = A;

if isfield(problem,'xhscale')
    xdc  = problem.x0scale(:)';
    xmax = problem.xhscale(:)';
    xscale = [xdc; repmat(xmax,hbm.harm.NFreq-1,1)*(1+1i)];
    problem.xscale = packdof(xscale,hbm.harm.iRetain);
    problem.Fscale = problem.xscale*0+1;
    problem.wscale = w;
    problem.Zscale = [problem.xscale; problem.wscale];
else
    problem = hbm_scaling(problem,hbm,init);
end

problem.Jscale = (1./problem.Fscale(:))*problem.Zscale(:)';

z0 = [x0; w];
Z0 = z0./problem.Zscale;

%and actually solve it
hbm.max_iter = 4;
bSuccess = false;
Jstr = [hbm.sparsity ones(NDof*NComp,1)];

constr_tol = 1E-6;
opt_tol = 1E-6;
maxit = 20;

attempts = 0;
while ~bSuccess && attempts < hbm.max_iter
    switch hbm.options.solver
        case 'fsolve'
            fun_obj = @(x)hbm_fsolve_obj(x,hbm,problem,A);
            fun_constr = @(x)hbm_fsolve_constr(x,hbm,problem,A);
            options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'Display','iter',...
                'OptimalityTolerance',opt_tol,'ConstraintTolerance',constr_tol,'MaxIterations',maxit);
            [Z,~,EXITFLAG,OUTPUT] = fmincon(fun_obj,Z0,[],[],[],[],[],[],fun_constr,options);
            bSuccess = EXITFLAG == 1;
            iter = OUTPUT.iterations + 1;
        case 'ipopt'
            options.jacobianstructure  = Jstr;
            options.jacobian = @hbm_jacobian;
            options.gradient = @hbm_grad;
            options.print_level = 5;
            options.max_iter = maxit;
            options.tol = opt_tol;
            options.constr_viol_tol = constr_tol;
            [Z, info] = fipopt(@hbm_obj,Z0,@hbm_constr,options,hbm,problem,A);           
            bSuccess = any(info.status == [0 1]);
            iter = info.iter;
    end
    Z0 = Z+rand(length(Z),1)*1E-8;
    attempts = attempts + 1;
end
if ~bSuccess
    Z = Z + NaN;
end
z = Z.*problem.Zscale;

w = abs(z(end));
w0 = w*hbm.harm.rFreqRatio;
x = z(1:end-1);

sol.w = w;
sol.A = A;
sol.X = unpackdof(x,hbm.harm.NFreq-1,NDof,hbm.harm.iRetain);
sol.U = A*feval(problem.excite,hbm,problem,w0);
sol.F = hbm_output3d(hbm,problem,w0,sol.U,sol.X);

%floquet multipliers & objective
u = packdof(sol.U);
sol.H = hbm_objective('complex',hbm,problem,w0,x,u);
sol.L = floquetMultipliers(hbm,problem,w0,u,x);

sol.it = iter;

function obj = hbm_obj(Z,hbm,problem,A)
x = Z(1:end-1).*problem.xscale;
w = Z(end).*problem.wscale;
w0 = w * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

H = hbm_objective('func',hbm,problem,w0,x,u);
obj = - problem.res.sign * H;

function G = hbm_grad(Z,hbm,problem,A)
x = Z(1:end-1).*problem.xscale;
w = Z(end).*problem.wscale;
w0 = w * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

[Dx, Dw] = hbm_objective({'jacobX','derivW'},hbm,problem,w0,x,u);
G = -problem.res.sign*[Dx Dw];

G = G.*problem.Zscale(:)';

function c = hbm_constr(Z,hbm,problem,A)
x = Z(1:end-1).*problem.xscale;
w = Z(end).*problem.wscale;
w0 = w * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

c = hbm_balance3d('func',hbm,problem,w0,u,x);

function J = hbm_jacobian(Z,hbm,problem,A)
x = Z(1:end-1).*problem.xscale;
w = Z(end).*problem.wscale;
w0 = w * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

Jx_nl = hbm_balance3d('jacob',hbm,problem,w0,u,x);
Dw_nl = hbm_balance3d('derivW',hbm,problem,w0,u,x);

J = [Jx_nl  Dw_nl];
J = J .* problem.Jscale;

function [obj,G] = hbm_fsolve_obj(Z,hbm,problem,A)
obj = hbm_obj(Z,hbm,problem,A);
if nargout > 1
    G =  hbm_grad(Z,hbm,problem,A)';
end

function [c,ceq,Jc,Jceq] = hbm_fsolve_constr(Z,hbm,problem,A)
ceq = hbm_constr(Z,hbm,problem,A);
c = [];
% c = ceq(end);
% ceq = ceq(1:end-1);
if nargout > 2
    Jceq = hbm_jacobian(Z,hbm,problem,A);
    Jc = zeros(0,length(Z));
    
    Jc = Jc';
    Jceq = Jceq';
%     Jc = Jceq(end,:)';
%     Jceq = Jceq(1:end-1,:)';
end

