function results = hbm_amp(hbm,problem,w0,A0,X0,AEnd,XEnd)
NDof = problem.NDof;

%first solve @ A0
sol = hbm_solve(hbm,problem,w0,A0,X0);
x0 = packdof(sol.X);
if any(isnan(abs(x0(:))))
    error('Failed to solve initial problem')
end
u0 = packdof(sol.U);

sol = hbm_solve(hbm,problem,w0,AEnd,XEnd);
xEnd = packdof(sol.X);
if any(isnan(abs(xEnd(:))))
    error('Failed to solve final problem')
end
uEnd = packdof(sol.U);

hbm.bIncludeNL = 1;
   
if isfield(problem,'xhscale')
    xdc  = problem.x0scale(:)';
    xmax = problem.xhscale(:)';
    xscale = [xdc; repmat(xmax,hbm.harm.NFreq-1,1)*(1+1i)];
    problem.xscale = packdof(xscale,hbm.harm.iRetain)*sqrt(length(xscale));
    problem.Ascale = mean([A0 AEnd]);
    problem.Fscale = problem.xscale*0+1;
    problem.Xscale = [problem.xscale; problem.Ascale];
    problem.Jscale = (1./problem.Fscale(:))*problem.Xscale(:)';
    bUpdateScaling = 0;
else
    problem = hbm_scaling(problem,hbm,x0,[],A0);
    bUpdateScaling = 1;
end
 
AMax = max(AEnd,A0);
AMin = min(AEnd,A0);

problem.AMax = AMax;
problem.AMin = AMin;
problem.A0   = A0;
problem.AEnd = AEnd;
problem.w0 = w0;

err = 'success';

switch hbm.cont.method
    case 'none'
        hbm_amp_plot('init',hbm,problem,x0,A0);
        A = A0;
        x = x0;
        u = u0;
        it = [];
        s = [];
        Ascale = mean([A0 AEnd]);
        
        hbm_amp_plot('data',hbm,problem,x0,A0);
        step = hbm.cont.step0;
        direction = sign(AEnd-A0)*Ascale;

        while A(end) <= AMax && A(end) >= AMin
            Apred = A(end) + step*direction;
            if Apred > AMax || Apred < AMin
                break;
            end
            
            iPredict = max(length(A)-6,1):length(A);
            if length(iPredict) > 1
                xpred = interp1(A(iPredict),x(:,iPredict)',Apred,'pchip','extrap')';
            else
                xpred = x;
            end
            
            %now try to solve
            Xpred = unpackdof(xpred,hbm.harm.NFreq-1,problem.NDof);
            sol = hbm_solve(hbm,problem,w0,Apred,Xpred);
            
            %unpack outputs           
            xsol = packdof(sol.X);
            usol = packdof(sol.U);
            
            if ~any(isnan(xsol))
                step = min(max(step * hbm.cont.C,hbm.cont.min_step),hbm.cont.max_step);
                x(:,end+1) = xsol;
                u(:,end+1) = usol;
                A(end+1) = Apred;
                it(end+1) = sol.it;
                s(end+1) = step;
                hbm_amp_plot('data',hbm,problem,xsol,Apred);
                num_err = 0;
                if A(end) >= AMax || A(end) <= AMin
                    break;
                end
            else
                step = step * hbm.cont.c;
                num_err = num_err + 1;
                hbm_amp_plot('err',hbm,problem,xsol,Apred);
                if num_err > hbm.cont.maxfail
                    err = 'Too many failed iterations';
                    break;
                end
            end
        end

        x(:,end+1) = xEnd;
        u(:,end+1) = uEnd;
        A(end+1)   = AEnd;
        it(end+1) = 0;
        s(end+1) = 0;
        
        hbm_amp_plot('close',hbm,problem,[],[]);
        
        debug = struct();
    case 'predcorr'
        x = x0; xCurr = x0;
        u = u0; uCurr = u0;
        A = A0; aCurr = A0;
        s = []; sCorrCurr = []; sPredCurr = [];
        it = []; itCurr = [];
        flag = {};
        
        step = hbm.cont.step0;
        num_fail = 0;
        num_step = 0;
        num_iter_tot = 0;
                        
        switch hbm.cont.predcorr.corrector
            case 'pseudo'
            case 'arclength'
                Jstr = [hbm.sparsity 0*x0+1;
                            0*x0'+1 1];
                switch hbm.cont.predcorr.solver
                    case 'ipopt'
                        ipopt_opt.print_level = 0;
                        ipopt_opt.maxit = hbm.cont.predcorr.maxit;
                        ipopt_opt.ftol = hbm.cont.ftol;
                        ipopt_opt.xtol = hbm.cont.xtol;
                        ipopt_opt.jacob = @hbm_arclength_jacobian;
                        ipopt_opt.jacobstructure = Jstr;
                    case 'fsolve'
                        fsolve_opt = optimoptions('fsolve',...
                            'Display','off',...
                            'FunctionTolerance',hbm.cont.ftol,...
                            'StepTolerance',hbm.cont.xtol,...
                            'SpecifyObjectiveGradient',true,...
                            'MaxIterations',hbm.cont.predcorr.maxit);
                end
        end
        
        hbm_amp_plot('init',hbm,problem,x,A);
        fprintf('STEP    PRED    CORR   STATUS  INFO     AMP     ITER   TOT      ')
        fprintf('X(%d)       ',1:length(x0))
        fprintf('\n')       
        
        fprintf('%3d   %6.4f  %6.4f    %s      %3s    %6.2f    %3d   %3d   ',0,0,0,'S','   ',A0,0,0)
        fprintf('%+5.2e  ',x0)
        fprintf('\n')
        
        Xprev = [x0; A0]./problem.Xscale;
        F = hbm_pseudo_constraints(Xprev,hbm,problem);
        J = hbm_pseudo_jacobian(Xprev,hbm,problem);
        tangent_prev = pseudo_null(J);
        tangent_prev = tangent_prev * sign(tangent_prev(end)) * sign(AEnd - A0);
        t0 = tangent_prev.*problem.Xscale;
        t = t0; tCurr = t0;
        Xend = [xEnd;AEnd]./problem.Xscale;
        
        while A(end) <= AMax && A(end) >= AMin
            
            %predictor
            switch hbm.cont.predcorr.predictor
                case 'linear'
                    Xpred = Xprev + step*tangent_prev;
                case 'quadratic'
                    if length(A) < 5
                        Xpred = Xprev + step*tangent_prev;
                    else
                    	Xpred = polynomial_predictor(Xsol(:,end-2:end),tsol(:,end-2:end),step);
                    end
                case 'cubic'
                    if length(A) < 5
                        Xpred = Xprev + step*tangent_prev;
                    else
                    	Xpred = polynomial_predictor(Xsol(:,end-3:end),tsol(:,end-3:end),step);
                    end
            end
            
            num_iter = 0;
            Xlast = Xpred + Inf;
            tangent = tangent_prev;
            
            %corrector
            switch hbm.cont.predcorr.corrector
                case 'pseudo'
                    bConverged = 0;
                    X = Xpred;
                    while  num_iter <= hbm.cont.predcorr.maxit
                        Xlast = X;
                        J = hbm_pseudo_jacobian(X,hbm,problem);
                        F = hbm_pseudo_constraints(X,hbm,problem);
                        if hbm.cont.predcorr.bMoorePenrose
                            X = Xlast - pinv(J)*F;
                            tangent = pseudo_null(J);
                        else
                            B = [J; tangent'];
                            R = [J*tangent; 0];
                            Q = [F; 0];
                            W = tangent - B\R;
                            tangent = W/norm(W);
                            X = Xlast - B\Q;
                        end
                        num_iter = num_iter + 1;
                       if ~(any(abs(X - Xlast) > hbm.cont.xtol) || any(abs(F) > hbm.cont.ftol))
                            bConverged = 1;
                            break;
                        end
                    end
                case 'arclength'
                    switch hbm.cont.predcorr.solver
                        case 'ipopt'
                            [X,info] = fipopt('',Xpred,@hbm_arclength_constraints,ipopt_opt,hbm,problem,Xprev,tangent_prev,step);
                            num_iter = info.iter;
                            bConverged = info.status == 0;
                            F = hbm_arclength_constraints(X,hbm,problem,Xprev,tangent_prev,step);
                        case 'fsolve'
                            [X,F,status,out] = fsolve(@hbm_arclength_constraints,Xpred,fsolve_opt,hbm,problem,Xprev,tangent_prev,step);
                            num_iter = out.iterations + 1;
                            bConverged = status == 1;
                    end
                    J = hbm_arclength_jacobian(X,hbm,problem,Xprev,tangent_prev,step);
                    tangent = pseudo_null(J(1:end-1,:));
                    tangent = sign(tangent'*tangent_prev)*tangent;
            end
            num_step = num_step + 1;
            num_iter_tot = num_iter_tot + num_iter;
            
            %prepare for plots
            aCurr(end+1) = X(end).*problem.Xscale(end);
            xCurr(:,end+1) = X(1:end-1).*problem.Xscale(1:end-1);
            uCurr(:,end+1) = packdof(aCurr(end)*feval(problem.excite,hbm,problem,w0*hbm.harm.rFreqRatio));
            sCorrCurr(end+1) = norm(X - Xprev)*sign((X-Xprev)'*tangent_prev);
            sPredCurr(end+1) = step;
            itCurr(end+1) = num_iter;
            
            step_prev = norm(X - Xprev);

            if bConverged && sign((X - Xprev)'*tangent_prev) > 0
                %success
                status = 'S';
                hbm_amp_plot('data',hbm,problem,xCurr(:,end),aCurr(end));
                                
                 if num_iter <= hbm.cont.num_iter_increase
                    step = min(sPredCurr(end) * hbm.cont.C,hbm.cont.max_step/sCorrCurr(end)*sPredCurr(end));
                    flag{end+1} = 'Success: Increasing step size';
                    info = 'Inc';
                elseif num_iter >= hbm.cont.num_iter_reduce
                    step = max(sPredCurr(end) / hbm.cont.C,hbm.cont.min_step/sCorrCurr(end)*sPredCurr(end));
                    flag{end+1} = 'Success: Reducing step size';
                    info = 'Red';
                else
                    flag{end+1} = 'Success';
                    info = '';
                end
                
                %just right
                
                if num_iter < hbm.cont.num_iter_reduce
                    step = min(step * hbm.cont.C,step * hbm.cont.max_step / step_prev);
                elseif num_iter > hbm.cont.num_iter_increase
                    step = step / hbm.cont.C;
                end

                %store the data
                A(end+1) = aCurr(end);
                x(:,end+1) = xCurr(:,end);
                u(:,end+1) = uCurr(:,end);
                t(:,end+1) = tCurr(:,end);
                s(end+1)   = sCorrCurr(end);
                it(end+1)  = itCurr(end);

                if bUpdateScaling
                    problem = hbm_scaling(problem,hbm,x(:,end),[],A(end));
                end

                Xsol = [x;A]./(repmat(problem.Xscale,1,size(x,2))); 
                tsol = normalise(t./(repmat(problem.Xscale,1,size(t,2))));
                
                Xend = [xEnd;AEnd]./problem.Xscale;
                
                Xprev = Xsol(:,end);
                tangent_prev = tsol(:,end);   
                
                num_fail = 0;
            else
                %failed
                hbm_amp_plot('err',hbm,problem,xCurr(:,end),aCurr(end));
                status = 'F';
                
                if w0 < 0
                    %gone to negative frequencies (somehow)
                    flag{end+1} = 'Failed: Negative frequency';
                    info = 'Neg';
                elseif num_iter <= hbm.cont.predcorr.maxit
                    %converge but step size is unacceptable
                    if step_prev < 0
                        %backwards
                        step = max(step * hbm.cont.c,hbm.cont.min_step);
                        flag{end+1} = 'Failed: Wrong direction';
                        info = 'Bwd';
                    elseif step_prev > hbm.cont.max_step
                        %too large
                        step = max(0.9 * step * hbm.cont.max_step / step_prev,hbm.cont.min_step);
                        flag{end+1} = 'Failed: Step too large';
                        info = 'Lrg';
                    elseif step_prev < hbm.cont.min_step
                        %too small
                        step = min(1.1 * step * hbm.cont.min_step / step_prev,hbm.cont.max_step);
                        flag{end+1} = 'Failed: Step too small';
                        info = 'Sml';
                    else
                        step = max(step * hbm.cont.c,hbm.cont.min_step);
                        flag{end+1} = 'Failed: Other error';
                        info = 'Oth';
                    end 
                else
                    %failed to converge
                    if any(abs(X - Xlast) > hbm.cont.xtol)
                        flag{end+1} = 'Failed: No convergence';
                        info = 'xtl';
                    elseif any(abs(F) > hbm.cont.ftol)
                        flag{end+1} = 'Failed: Constraints violated';
                        info = 'ftl';
                    else
                        flag{end+1} = 'Failed';
                        info = '';
                    end
                    step = max(step * hbm.cont.c,hbm.cont.min_step);
                end
                
                num_fail = num_fail + 1;
                if num_fail > hbm.cont.maxfail
                    err = 'Too many failed iterations';
                    break
                elseif w0 < 0
                    err = 'Negative frequency';
                    break
                end
            end
            fprintf('%3d   %6.4f  %6.4f    %s      %3s    %6.2f    %3d   %3d   ',num_step,step,step_prev,status,info,aCurr(end),num_iter,num_iter_tot)
            fprintf('%+5.2e  ',xCurr(:,end))
            fprintf('\n')
        end
        
        %add on final point
        aCurr(end+1)   = AEnd;
        xCurr(:,end+1) = xEnd;
        uCurr(:,end+1) = uEnd;
        sCorrCurr(end+1) = norm(Xend - Xprev);
        sPredCurr(end+1) = norm(Xend - Xprev);
        flag{end+1} = 'Success';
        itCurr(end+1) = 0;
        
        %store the data
        A(end+1)   = aCurr(end);
        x(:,end+1) = xCurr(:,end);
        u(:,end+1) = uCurr(:,end);
        s(end+1)   = sCorrCurr(end);
        it(end+1)  = itCurr(end);
        
        hbm_amp_plot('close',hbm,problem,[],[]);
        debug = struct('x',xCurr,...
            'a',aCurr,...
            'u',uCurr,...
            'sCorr',sCorrCurr,...
            'sPred',sPredCurr,...
            'it',itCurr,...
            'flag',flag);
    case 'coco'
        rng('shuffle')
        currdir = pwd;
        cd(desktoproot);
        name = ['temp' num2str(round(rand(1)*1E12))];
        while exist(['data' filesep name],'dir')
            name = ['temp' num2str(round(rand(1)*1E12))];
        end
        data.hbm = hbm;
        data.problem = problem;
                
        prob = coco_prob();
        prob = coco_add_func(prob, 'alg', @hbm_coco_constraints, @hbm_coco_jacobian, data, 'zero','u0', [x0(:); A0]./problem.Xscale);
        prob = coco_add_pars(prob, 'pars', length(x0)+1, 'A');
        prob = coco_add_slot(prob, 'plot', @hbm_coco_callback, data, 'bddat');

        prob = coco_set(prob,'cont','h0',hbm.cont.step0,'h_min',hbm.cont.min_step,'h_max',hbm.cont.max_step,'ItMX',hbm.cont.coco.ItMX,'NPR',hbm.cont.coco.NPR);
        prob = coco_set(prob,'ep','bifus',false);
        bd = coco(prob, 'temp', [], 1, 'A', [AMin AMax]./problem.Xscale(end));
        hbm_amp_plot('close',hbm,problem,[],[]);
        
        %extract the solutions
        lab_col = coco_bd_col(bd, 'TYPE');
        NPts = sum(cellfun(@(x)~isempty(x),lab_col));
        x = zeros(hbm.harm.NComp*NDof,NPts);
        A = zeros(1,NPts);
        for i = 1:NPts
            chart = coco_read_solution('temp',i,'chart');
            x(:,i) = chart.x(1:end-2);
            A(i) = chart.x(end);
        end
        A = A .* problem.Ascale;
        x = x .* (problem.xscale*(0*A+1));
        fclose all;
        try
            rmdir('data','s')
        end
        cd(currdir);
        %TODO: extract the paramaters from coco
        debug = struct();
        s = 0*A;
        it = 0*A;
        err = 'success';
end

X = unpackdof(x,hbm.harm.NFreq-1,NDof);
w = w0 + 0*A;

W = (hbm.harm.kHarm*(hbm.harm.rFreqRatio.*hbm.harm.rFreqBase)')*w;

NPts = size(x,2);
if ~exist('u','var')
    u = zeros(hbm.harm.NComp*problem.NInput,NPts);
    for i = 1:NPts
        w0 = w(i)*hbm.harm.rFreqRatio;
        U = A(i)*feval(problem.excite,hbm,problem,w0);
        u(:,i) = packdof(U);
    end
end
U = unpackdof(u,hbm.harm.NFreq-1,problem.NInput);

if ~exist('f','var')
    f = zeros(hbm.harm.NComp*problem.NOutput,NPts);
    for i = 1:NPts
        w0 = w(i)*hbm.harm.rFreqRatio;
        f(:,i) = hbm_output3d(hbm,problem,w0,u(:,i),x(:,i));
    end
end
F = unpackdof(f,hbm.harm.NFreq-1,problem.NOutput);

for i = 1:NPts
    w0 = w(i)*hbm.harm.rFreqRatio;
    lambda(:,i) = floquetMultipliers(hbm,problem,w0,u(:,i),x(:,i));
end

results = struct('X',X,...
    'U',U,...
    'F',F,...
    'w',w,...
    'W',W,...
    'A',A,...
    's',s,...
    'it',it,...
    'L',lambda,...
    'err',err,...
    'debug',debug);

%% Pseudo functions
function y = norm2(x)
y = sqrt(sum(x.^2,1));

function X_extrap = polynomial_predictor(X,dX,s_extrap)
s = norm2(diff(X,[],2));
s = cumsum([0 s]);
N = length(s);
if ~isempty(dX)
    [A,Ad] = poly_mat(s,N);
    p = [A;Ad]\[X dX]';
else
    A = poly_mat(s,N);
    p = A\X';
end
B = poly_mat(s(end) + s_extrap,N);
X_extrap = (B*p)';

function [A,Ad] = poly_mat(s,N)
A = zeros(length(s),N);
Ad = zeros(length(s),N);
for i = 1:N
    A(:,i) = s.^(i-1);
    if nargout > 1
        if i > 1
            Ad(:,i) = (i-1)*s.^(i-2);
        else
            Ad(:,i) = 0*s;
        end
    end
end

function t = pseudo_null(A)
[U,S,V] = svd(A,0);
t = V(:,end);
t = t./norm(t);

function c = hbm_pseudo_constraints(X,hbm,problem)
%unpack the inputs
Xscale = problem.Xscale;
x = X(1:end-1).*Xscale(1:end-1);
A = X(end).*Xscale(end);
w0 = problem.w0 * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

c = hbm_balance3d('func',hbm,problem,w0,u,x);

function J = hbm_pseudo_jacobian(X,hbm,problem)

%unpack the inputs
Xscale = problem.Xscale;
x = X(1:end-1).*Xscale(1:end-1);
A = X(end).*Xscale(end);
w0 = problem.w0 * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

Jx = hbm_balance3d('jacob',hbm,problem,w0,u,x);
Da = hbm_balance3d('derivA',hbm,problem,w0,u,x);

J = [Jx  Da];

J = J .* repmat(Xscale(:)',size(J,1),1);


%% Coco functions

function [data res] = hbm_coco_callback(prob, data, command, varargin)
hbm = data.hbm;
problem = data.problem;

switch command
    case 'init'
        res = {};
        x0 = prob.efunc.x0(1:end-1).*problem.xscale;
        A0 = prob.efunc.x0(end).*problem.Ascale;
        hbm_amp_plot('init',hbm,problem,x0,A0);
    case 'data'
        chart = varargin{1};
        x = chart.x(1:end-2).*problem.xscale;
        a = chart.x(end).*problem.Ascale;
        hbm_amp_plot('data',hbm,problem,x,a);
        res = {};
end

function [data,c] = hbm_coco_constraints(prob, data, u)
hbm = data.hbm;
problem = data.problem;
Xscale = problem.Xscale;

x = u(1:end-1).*Xscale(1:end-1);
A = u(end).*Xscale(end);
w0 = problem.w0 * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

c = hbm_balance3d('func',hbm,problem,w0,u,x);

function [data, J] = hbm_coco_jacobian(prob, data, u)
hbm = data.hbm;
problem = data.problem;
Xscale = problem.Xscale;

x = u(1:end-1).*Xscale(1:end-1);
A = u(end).*Xscale(end);
w0 = problem.w0 * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

Jx = hbm_balance3d('jacob',hbm,problem,w0,u,x);
Da = hbm_balance3d('derivA',hbm,problem,w0,u,x)/A;

J = [Jx  Da];

J = J .* repmat(Xscale(:)',size(J,1),1);

%% Arclength files

function [f,J] = hbm_arclength_constraints(X,hbm,problem,Xprev,tangent_prev,step)
Xscale = problem.Xscale;
Xcurr = X;

x = X(1:end-1).*Xscale(1:end-1);
A = X(end).*Xscale(end);
w0 = problem.w0 * hbm.harm.rFreqRatio;

X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof);
U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

c = hbm_balance3d('func',hbm,problem,w0,u,x);

sgn = sign((Xcurr - Xprev)' * tangent_prev);
s = norm(Xcurr - Xprev) * sgn;

f = [c;
     s - step];
 
if nargout > 1
    J = hbm_arclength_jacobian(Xcurr,hbm,problem,Xprev,tangent_prev,step);
end
      
function J = hbm_arclength_jacobian(X,hbm,problem,Xprev,tangent_prev,step)
Xscale = problem.Xscale;
Xcurr = X;

x = X(1:end-1).*Xscale(1:end-1);
A = X(end).*Xscale(end);
w0 = problem.w0 * hbm.harm.rFreqRatio;

X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof);
U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

sgn = sign((Xcurr - Xprev)' * tangent_prev);

Jx = hbm_balance3d('jacob',hbm,problem,w0,u,x);
Da = hbm_balance3d('derivA',hbm,problem,w0,u,x)/A;

J = [Jx  Da];

J = J .* repmat(Xscale(:)',size(J,1),1);

J(end+1,:) = sgn*((Xcurr - Xprev)'+eps)./(norm(Xcurr - Xprev)+eps);