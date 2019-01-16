function results = hbm_frf(hbm,problem,A,w0,x0,wEnd,xEnd)
NDof = problem.NDof;

%first solve @ w0
sol = hbm_solve(hbm,problem,w0,A,x0);
x0 = packdof(sol.X,hbm.harm.iRetain);
if any(isnan(abs(x0(:))))
    error('Failed to solve initial problem')
end
u0 = packdof(sol.U);
f0 = packdof(sol.F);

sol = hbm_solve(hbm,problem,wEnd,A,xEnd);
xEnd = packdof(sol.X,hbm.harm.iRetain);
if any(isnan(abs(xEnd(:))))
    error('Failed to solve final problem')
end
uEnd = packdof(sol.U);
fEnd = packdof(sol.F);

hbm.bIncludeNL = 1;

if isfield(problem,'xhscale')
    xdc  = problem.x0scale(:)';
    xmax = problem.xhscale(:)';
    xscale = [xdc; repmat(xmax,hbm.harm.NFreq-1,1)*(1+1i)];
    problem.xscale = packdof(xscale,hbm.harm.iRetain)*sqrt(length(xscale));
    problem.wscale = mean([w0 wEnd]);
    problem.Fscale = problem.xscale*0+1;
    problem.Xscale = [problem.xscale; problem.wscale];
    problem.Jscale = (1./problem.Fscale(:))*problem.Xscale(:)';
    bUpdateScaling = 0;
else
    problem = hbm_scaling(problem,hbm,x0,w0);
    bUpdateScaling = 1;
end

wMax = max(wEnd,w0);
wMin = min(wEnd,w0);

problem.wMax = wMax;
problem.wMin = wMin;
problem.w0   = w0;
problem.wEnd = wEnd;

err = 'success';

switch hbm.cont.method
    case 'none'
        hbm_frf_plot('init',hbm,problem,x0,w0,A);
        w = w0;
        x = x0;
        u = u0;
        it = [];
        s = [];
        wscale = mean([w0 wEnd]);
        
        hbm_frf_plot('data',hbm,problem,x0,w0,A);
        step = hbm.cont.step0;
        direction = sign(wEnd-w0)*wscale;
        
        while  w(end) >= wMin && w(end) <= wMax
            wpred = w(end) + step*direction;
            if wpred > wMax || wpred < wMin
                break;
            end
            
            iPredict = max(length(w)-6,1):length(w);
            if length(iPredict) > 1
                xpred = interp1(w(iPredict),x(:,iPredict)',wpred,'pchip','extrap')';
            else
                xpred = x;
            end
 
            Xpred = unpackdof(xpred,hbm.harm.NFreq-1,problem.NDof,hbm.harm.iRetain);
            sol = hbm_solve(hbm,problem,wpred,A,Xpred);
            
            %unpack outputs           
            xsol = packdof(sol.X);
            usol = packdof(sol.U);
            wsol = sol.w0;
            
            if ~any(isnan(xsol))
                step = min(max(step * hbm.cont.C,hbm.cont.min_step),hbm.cont.max_step);
                x(:,end+1) = xsol;
                u(:,end+1) = usol;
                w(end+1) = wsol;
                it(end+1) = sol.it;
                s(end+1) = step;
                hbm_frf_plot('data',hbm,problem,xsol,wsol,A);
                num_err = 0;
                if w(end) >= wMax || w(end) <= wMin
                    break;
                end
            else
                step = step * hbm.cont.c;
                num_err = num_err + 1;
                hbm_frf_plot('err',hbm,problem,xsol,wsol,A);
                if num_err > hbm.cont.maxfail
                    err = 'Too many failed iterations';
                    break;
                end
            end
        end
        
        x(:,end+1) = xEnd;
        u(:,end+1) = uEnd;
        w(end+1)   = wEnd;
        it(end+1)  = 0;
        s(end+1) = 0;
        
        hbm_frf_plot('close',hbm,problem,[],[],[]);
        
        debug = struct();
    case 'predcorr'              
        w = w0; wCurr = w0;
        x = x0; xCurr = x0;
        u = u0; uCurr = u0;
        f = f0; fCurr = f0;
        s = []; sCorrCurr = []; sPredCurr = [];
        it = []; itCurr = [];
        flag = {};
        
        Xprev = [x0; w0]./problem.Xscale;

        num_step = 1;
        num_iter = 0;
        num_fail = 0;
        num_iter_tot = 0;
        
        step = hbm.cont.step0;
        
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
        
        hbm_frf_plot('init',hbm,problem,x,w,A);
        fprintf('STEP    PRED    CORR   STATUS  INFO  ITER   TOT   FREQ       ')
        fprintf('X(%d)       ',1:length(x0))
        fprintf('\n')
        fprintf('%3d   %6.4f  %6.4f    %s       %3s  %3d   %3d   %6.2f   ',num_step,step,step,'S','   ',num_iter,num_iter_tot,w)
        fprintf('%+5.2e  ',x)
        fprintf('\n')       
        
        Xprev = [x0; w0]./problem.Xscale;
        F = hbm_pseudo_constraints(Xprev,hbm,problem,A);
        J = hbm_pseudo_jacobian(Xprev,hbm,problem,A);
        tangent_prev = pseudo_null(J);
        tangent_prev = tangent_prev * sign(tangent_prev(end)) * sign(wEnd - w0);
        t0 = tangent_prev.*problem.Xscale;
        t = t0; tCurr = t0;
        Xend = [xEnd;wEnd]./problem.Xscale;
        
        while norm(Xprev - Xend) > hbm.cont.max_step
            
            %predictor
            switch hbm.cont.predcorr.predictor
                case 'linear'
                    Xpred = Xprev + step*tangent_prev;
                case 'quadratic'
                    if length(w) < 5
                        Xpred = Xprev + step*tangent_prev;
                    else
                    	Xpred = polynomial_predictor(Xsol(:,end-2:end),tsol(:,end-2:end),step);
                    end
                case 'cubic'
                    if length(w) < 5
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
                    while num_iter <= hbm.cont.predcorr.maxit
                        Xlast = X;
                        J = hbm_pseudo_jacobian(X,hbm,problem,A);
                        F = hbm_pseudo_constraints(X,hbm,problem,A);
                        if hbm.cont.predcorr.bMoorePenrose
                            X = Xlast - J\F;
                        else
                            B = [J; tangent'];
                            R = [J*tangent; 0];
                            Q = [F; 0];
                            W = tangent - B\R;
                            tangent = normalise(W);
                            X = Xlast - B\Q;
                        end
                        num_iter = num_iter + 1;
                        if ~(any(abs(X - Xlast) > hbm.cont.xtol) || any(abs(F) > hbm.cont.ftol))
                            bConverged = 1;
                            break;
                        end
                    end
                    if hbm.cont.predcorr.bMoorePenrose
                        tangent = pseudo_null(J);
                        tangent = sign(tangent'*tangent_prev)*tangent;
                    end
                case 'arclength'
                    switch hbm.cont.predcorr.solver
                        case 'ipopt'
                            [X,info] = fipopt('',Xpred,@hbm_arclength_constraints,ipopt_opt,hbm,problem,A,Xprev,tangent_prev,step);
                            num_iter = info.iter;
                            bConverged = info.status == 0;
                            F = hbm_arclength_constraints(X,hbm,problem,A,Xprev,tangent_prev,step);
                        case 'fsolve'
                            [X,F,status,out] = fsolve(@hbm_arclength_constraints,Xpred,fsolve_opt,hbm,problem,A,Xprev,tangent_prev,step);
                            num_iter = out.iterations + 1;
                            bConverged = status == 1;
                    end
                    J = hbm_arclength_jacobian(X,hbm,problem,A,Xprev,tangent_prev,step);
                    tangent = pseudo_null(J(1:end-1,:));
                    tangent = sign(tangent'*tangent_prev)*tangent;
            end

            num_step = num_step + 1;
            num_iter_tot = num_iter_tot + num_iter;

            %prepare for plots
            wCurr(end+1) = X(end).*problem.wscale;
            xCurr(:,end+1) = X(1:end-1).*problem.xscale;
            tCurr(:,end+1) = tangent.*problem.Xscale;
            uCurr(:,end+1) = packdof(A*feval(problem.excite,hbm,problem,wCurr(end)*hbm.harm.rFreqRatio));
            fCurr(:,end+1) = hbm_output3d(hbm,problem,wCurr(end)*hbm.harm.rFreqRatio,uCurr(:,end),xCurr(:,end));
            sCorrCurr(end+1) = norm(X - Xprev)*sign((X-Xprev)'*tangent_prev);
            sPredCurr(end+1) = step;
            itCurr(end+1) = num_iter;
            
            if bConverged && sCorrCurr(end) >= hbm.cont.min_step && sCorrCurr(end) <= hbm.cont.max_step && wCurr(end) > 0
                %sucess
                status = 'S';
                hbm_frf_plot('data',hbm,problem,xCurr(:,end),wCurr(end),A);
                             
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
                
                %store the data
                w(end+1)   = wCurr(end);
                x(:,end+1) = xCurr(:,end);
                u(:,end+1) = uCurr(:,end);
                t(:,end+1) = tCurr(:,end);
                f(:,end+1) = fCurr(:,end);
                s(end+1)   = sCorrCurr(end);
                it(end+1)  = itCurr(end);
                
                if bUpdateScaling
                    problem = hbm_scaling(problem,hbm,x(:,end),w(end));
                end

                Xsol = [x;w]./(repmat(problem.Xscale,1,size(x,2))); 
                tsol = normalise(t./(repmat(problem.Xscale,1,size(t,2))));

                Xend = [xEnd;wEnd]./problem.Xscale;

                Xprev = Xsol(:,end);
                tangent_prev = tsol(:,end);              

                num_fail = 0;
            else
                hbm_frf_plot('err',hbm,problem,xCurr(:,end),wCurr(end),A);
                status = 'F';
                
                if wCurr(end) < 0
                    %gone to negative frequencies (somehow)
                    flag{end+1} = 'Failed: Negative frequency';
                    info = 'Neg';
                elseif num_iter <= hbm.cont.predcorr.maxit
                    %converge but step size is unacceptable
                    if sCorrCurr(end) < 0
                        %backwards
                        step = max(sPredCurr(end) * hbm.cont.c,hbm.cont.min_step);
                        flag{end+1} = 'Failed: Wrong direction';
                        info = 'Bwd';
                    elseif sCorrCurr(end) > hbm.cont.max_step
                        %too large
                        step = max(0.9 * sPredCurr(end) * hbm.cont.max_step / sCorrCurr(end),hbm.cont.min_step);
                        flag{end+1} = 'Failed: Step too large';
                        info = 'Lrg';
                    elseif sCorrCurr(end) < hbm.cont.min_step
                        %too small
                        step = min(1.1 * sPredCurr(end) * hbm.cont.min_step / sCorrCurr(end),hbm.cont.max_step);
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
                    step = max(sPredCurr(end) * hbm.cont.c,hbm.cont.min_step);
                end
                
                num_fail = num_fail + 1;
                if num_fail > hbm.cont.maxfail
                    err = 'Too many failed iterations';
                    break
                elseif wCurr(end) < 0
                    err = 'Negative frequency';
                    break
                end
            end
            
            fprintf('%3d   %6.4f  %6.4f    %s       %3s  %3d   %3d   %6.2f   ',num_step,sPredCurr(end),sCorrCurr(end),status,info,num_iter,num_iter_tot,wCurr(end))
            fprintf('%+5.2e  ',xCurr(:,end))
            fprintf('\n')
        end
        
        %add on final point
        wCurr(end+1)   = wEnd;
        xCurr(:,end+1) = xEnd;
        uCurr(:,end+1) = uEnd;
        fCurr(:,end+1) = fEnd;
        sCorrCurr(end+1) = norm(Xend - Xprev);
        sPredCurr(end+1) = norm(Xend - Xprev);
        flag{end+1} = 'Success';
        itCurr(end+1) = 0;
        
        %store the data
        w(end+1)   = wCurr(end);
        x(:,end+1) = xCurr(:,end);
        u(:,end+1) = uCurr(:,end);
        f(:,end+1) = fCurr(:,end);
        t(:,end+1) = tCurr(:,end);
        s(end+1)   = sCorrCurr(end);
        it(end+1)  = itCurr(end);

        hbm_frf_plot('close',hbm,problem,[],[],[]);
        debug = struct('x',xCurr,...
            'u',uCurr,...
            'f',fCurr,...
            'w',wCurr,...
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
        data.A = A;

        prob = coco_prob();
        prob = coco_add_func(prob, 'alg', @hbm_coco_constraints, @hbm_coco_jacobian, data, 'zero','u0', [x0(:); w0]./problem.Xscale);
        prob = coco_add_pars(prob, 'pars', length(x0)+1, 'w');
        prob = coco_add_slot(prob, 'plot', @hbm_coco_callback, data, 'bddat');
        
        prob = coco_set(prob,'cont','h0',hbm.cont.step0,'h_min',hbm.cont.min_step,'h_max',hbm.cont.max_step,'ItMX',hbm.cont.coco.ItMX,'NPR',hbm.cont.coco.NPR);
        prob = coco_set(prob,'ep','bifus',false);
        bd = coco(prob, name, [], 1, 'w', [wMin wMax]./problem.wscale);
        hbm_frf_plot('close',hbm,problem,[],[],[]);
        
        %extract the solutions
        lab_col = coco_bd_col(bd, 'TYPE');
        NPts = sum(cellfun(@(x)~isempty(x),lab_col));
        x = zeros(hbm.harm.NComp*NDof,NPts);
        w = zeros(1,NPts);
        for i = 1:NPts
            chart = coco_read_solution(name,i,'chart');
            x(:,i) = chart.x(1:end-2);
            w(i) = chart.x(end);
        end
        w = w .* problem.wscale;
        x = x .* (problem.xscale*(0*w+1));
        fclose all;
        try
            rmdir(['data' filesep name],'s')
        end
        cd(currdir);
        %TODO: extract the paramaters from coco
        debug = struct();
        s = 0*w;
        it = 0*w;
        err = 'success';
end

X = unpackdof(x,hbm.harm.NFreq-1,NDof,hbm.harm.iRetain);

W = (hbm.harm.kHarm*(hbm.harm.rFreqRatio.*hbm.harm.rFreqBase)')*w;

NPts = size(x,2);
if ~exist('u','var')
    u = zeros(hbm.harm.NComp*problem.NInput,NPts);
    for i = 1:NPts
        w0 = w(i)*hbm.harm.rFreqRatio;
        U = A*feval(problem.excite,hbm,problem,w0);
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

if ~exist('lambda','var')
    lambda = zeros(2*size(x,1),NPts);
    for i = 1:NPts
        w0 = w(i)*hbm.harm.rFreqRatio;
        lambda(:,i) = floquetMultipliers(hbm,problem,w0,u(:,i),x(:,i));
    end
end

A = A+0*w;
results = struct('X',X,...
    'A',A,...
    'U',U,...
    'F',F,...
    'w',w,...
    'W',W,...
    's',s,...
    'it',it,...
    'L',lambda,...
    'err',err,...
    'debug',debug);

function y = norm2(x)
y = sqrt(sum(x.^2,1));

function y = normalise(x)
scale = repmat(norm2(x),size(x,1),1);
y = x ./scale;

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

%% Pseudo functions

function t = pseudo_null(A)
[U,S,V] = svd(A,0);
t = V(:,end);
t = t./norm(t);

function c = hbm_pseudo_constraints(X,hbm,problem,A)
%unpack the inputs
x = X(1:end-1).*problem.xscale;
w0 = X(end).*problem.wscale;

w0 = w0 * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

c = hbm_balance3d('func',hbm,problem,w0,u,x);
c = c ./ problem.Fscale;

function J = hbm_pseudo_jacobian(X,hbm,problem,A)
% J = jacobian(@hbm_pseudo_constraints,X,hbm,problem,A);
% return;

%unpack the inputs
x = X(1:end-1).*problem.xscale;
w0 = X(end).*problem.wscale;

w0 = w0 * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

Jx = hbm_balance3d('jacob',hbm,problem,w0,u,x);
Dw = hbm_balance3d('derivW',hbm,problem,w0,u,x);

J = [Jx Dw];

% hbm.options.bAnalyticalDerivs = false;
% Jx2 = hbm_balance3d('jacob',hbm,problem,w0,u,x);
% Dw2 = hbm_balance3d('derivW',hbm,problem,w0,u,x);
% J2 = [Jx2 Dw2];

J = J .* problem.Jscale;


%% Arclength files
function [f,J] = hbm_arclength_constraints(X,hbm,problem,A,Xprev,tangent_prev,step)
x = X(1:end-1).*problem.xscale;
w0 = X(end).*problem.wscale;

w0 = w0 * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

c = hbm_balance3d('func',hbm,problem,w0,u,x);
c = c ./ problem.Fscale;

sgn = sign((X - Xprev)' * tangent_prev);
s = norm(X - Xprev) * sgn;
f = [c;
    s - step];

if nargout > 1
    J = hbm_arclength_jacobian(X,hbm,problem,A,Xprev,tangent_prev,step);
end

function J = hbm_arclength_jacobian(X,hbm,problem,A,Xprev,tangent_prev,step)
x = X(1:end-1).*problem.xscale;
w0 = X(end).*problem.wscale;

w0 = w0 * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

sgn = sign((X - Xprev)' * tangent_prev);

Jx = hbm_balance3d('jacob',hbm,problem,w0,u,x);
Dw = hbm_balance3d('derivW',hbm,problem,w0,u,x);

J = [Jx Dw];
J = J .*problem.Jscale;

J(end+1,:) = sgn*((X - Xprev)'+eps)/(1*(norm(X - Xprev)+eps)); %length(X)

%% Coco functions

function [data res] = hbm_coco_callback(prob, data, command, varargin)
hbm = data.hbm;
problem = data.problem;
A = data.A;
step = prob.cont.arc_alpha;

switch command
    case 'init'
        res = {};
        x0 = prob.efunc.x0(1:end-1).*problem.xscale;
        w0 = prob.efunc.x0(end).*problem.wscale;
        hbm_frf_plot('init',hbm,problem,x0,w0,A);
    case 'data'
        chart = varargin{1};
        x = chart.x(1:end-2).*problem.xscale;
        w = chart.x(end).*problem.wscale;
        hbm_frf_plot('data',hbm,problem,x,w,A);
        res = {};
end

function [data,c] = hbm_coco_constraints(prob, data, u)
hbm = data.hbm;
problem = data.problem;
A = data.A;

x = u(1:end-1).*problem.xscale;
w0 = u(end).*problem.wscale;

w0 = w0 * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

c = hbm_balance3d('func',hbm,problem,w0,u,x);
c = c ./ problem.Fscale;


function [data, J] = hbm_coco_jacobian(prob, data, u)
hbm = data.hbm;
problem = data.problem;
A = data.A;

x = u(1:end-1).*problem.xscale;
w0 = u(end).*problem.wscale;

w0 = w0 * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

%multiple fundamental
Jx = hbm_balance3d('jacob',hbm,problem,w0,u,x);
Dw = hbm_balance3d('derivW',hbm,problem,w0,u,x);

J = [Jx Dw];

% hbm.options.bAnalyticalDerivs = false;
% Jx2 = hbm_balance('jacob',hbm,problem,w0,u,x);
% Dw2 = hbm_balance('derivW',hbm,problem,w0,u,x);
% J2 = [Jx2 Dw2];
%
% if max(abs(Dw -Dw2)) > 1E-4
%     keyboard
% end
%
% if max(abs(Jx(:) - Jx2(:))) > 1E-4
%     keyboard
% end

J = J .*problem.Jscale;