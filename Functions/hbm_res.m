function sol = hbm_res(hbm,problem,w,A,X0)
problem.type = 'resonance';

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
hbm.bIncludeNL = true;

%first solve @ w0
sol = hbm_solve(hbm,problem,w,A,X0);
X0 = sol.X;
z0 = packdof(X0(:,problem.iNL));

init.X = X0;
init.w = w;
init.A = A;

if isfield(problem,'xscale')
    xscale = [problem.xscale'; repmat(problem.xscale',hbm.harm.NFreq-1,1)*(1+1i)];   
    problem.Xscale = packdof(xscale(:,problem.iNL))*sqrt(length(xscale));
    problem.wscale = w;
    problem.Fscale = problem.Xscale*0+1;
    problem.Zscale = [problem.Xscale; problem.wscale];
    problem.Jscale = (1./problem.Fscale(:))*problem.Zscale(:)';
else
    problem = hbm_scaling(problem,hbm,init);
end

%and actually solve it
hbm.max_iter = 4;
bSuccess = false;

constr_tol = 1E-6;
opt_tol = 1E-6;
maxit = 20;

Z0 = [z0; w]./problem.Zscale;

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
            options.jacobianstructure  = [hbm.sparsity ones(problem.NNL*hbm.harm.NComp,1)];
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

w = z(end);
w0 = w*hbm.harm.rFreqRatio + hbm.harm.wFreq0;
W = hbm.harm.kHarm*(hbm.harm.rFreqBase.*w0)';

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

x = hbm_recover(hbm,problem,w,u,z(1:end-1));
X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof);

sol.w = w;
sol.W = W;
sol.A = A;
sol.X = X;
sol.U = U;
sol.F = hbm_output3d(hbm,problem,w,sol.U,sol.X);

%floquet multipliers & objective
sol.H = hbm_objective('complex',hbm,problem,w,z(1:end-1),u);
sol.L = hbm_floquet(hbm,problem,w,u,x);

sol.it = iter;

function obj = hbm_obj(Z,hbm,problem,A)
x = Z(1:end-1).*problem.Xscale;
w = Z(end).*problem.wscale;
w0 = w * hbm.harm.rFreqRatio + hbm.harm.wFreq0;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

H = hbm_objective('func',hbm,problem,w,x,u);
obj = - problem.res.sign * H;

function G = hbm_grad(Z,hbm,problem,A)
x = Z(1:end-1).*problem.Xscale;
w = Z(end).*problem.wscale;
w0 = w * hbm.harm.rFreqRatio + hbm.harm.wFreq0;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

[Dx, Dw] = hbm_objective({'jacobX','derivW'},hbm,problem,w,x,u);
G = -problem.res.sign*[Dx Dw];

G = G.*problem.Zscale(:)';

function c = hbm_constr(Z,hbm,problem,A)
x = Z(1:end-1).*problem.Xscale;
w = Z(end).*problem.wscale;
w0 = w * hbm.harm.rFreqRatio + hbm.harm.wFreq0;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

c = hbm_balance3d('func',hbm,problem,w,u,x);

function J = hbm_jacobian(Z,hbm,problem,A)
x = Z(1:end-1).*problem.Xscale;
w = Z(end).*problem.wscale;
w0 = w * hbm.harm.rFreqRatio + hbm.harm.wFreq0;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

Jx_nl = hbm_balance3d('jacob',hbm,problem,w,u,x);
Dw_nl = hbm_balance3d('derivW',hbm,problem,w,u,x);

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