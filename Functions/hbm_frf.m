function [results,curr] = hbm_frf(hbm,problem,A,w0,X0,wEnd,XEnd)
problem.type = 'frf';
problem.A = A;

%first solve @ w0
sol = hbm_solve(hbm,problem,w0,A,X0);
x0 = packdof(sol.X(:,problem.iNL),hbm.harm.iRetainNL);
u0 = packdof(sol.U);
f0 = packdof(sol.F);
if any(isnan(abs(x0(:))))
    results = struct('X',x0,...
    'A',A,...
    'U',u0,...
    'F',f0,...
    'w',w0,...
    'err','Failed to solve initial problem');
    return;
end
z0 = [x0;w0];

init.X = sol.X;
init.w = w0;
init.A = A;

sol = hbm_solve(hbm,problem,wEnd,A,XEnd);
xEnd = packdof(sol.X(:,problem.iNL),hbm.harm.iRetainNL);
uEnd = packdof(sol.U);
fEnd = packdof(sol.F);
if any(isnan(abs(xEnd(:))))
    results = struct('X',xEnd,...
    'A',A,...
    'U',uEnd,...
    'F',fEnd,...
    'w',wEnd,...
    'err','Failed to solve final problem');
    return;
end
zEnd = [xEnd;wEnd];

hbm.bIncludeNL = 1;

if isfield(problem,'xhscale')
    xdc  = problem.x0scale(:)';
    xmax = problem.xhscale(:)';
    xscale = [xdc; repmat(xmax,hbm.harm.NFreq-1,1)*(1+1i)];
    problem.xscale = packdof(xscale,hbm.harm.iRetain)*sqrt(length(xscale));
    problem.wscale = mean([w0 wEnd]);
    problem.Fscale = problem.xscale*0+1;
    problem.Zscale = [problem.xscale; problem.wscale];
    problem.Jscale = (1./problem.Fscale(:))*problem.Zscale(:)';
    bUpdateScaling = 0;
else
    problem = hbm_scaling(problem,hbm,init);
    bUpdateScaling = 1;
end

wMax = max(wEnd,w0);
wMin = min(wEnd,w0);

problem.wMax = wMax;
problem.wMin = wMin;
problem.w0   = w0;
problem.wEnd = wEnd;

prog.Status = 'success';

switch hbm.cont.method
    case 'none'
        hbm_frf_plot('init',hbm,problem,init);
        wscale = mean([w0 wEnd]);
        
        pred.step = hbm.cont.step0;
        direction = sign(wEnd-w0)*wscale;
        
        problem.xscale = 0*problem.xscale + 1;
        problem.wscale = 1;
        problem.Zscale = 0*problem.Zscale + 1;
        
        z = z0;
        J = hbm_frf_jacobian(z,hbm,problem);
        t = get_tangent(J);
        
        corr.step = 0;
        corr.it = 0;
        
        curr = hbm_frf_results(z,t,pred,corr,hbm,problem);
        results = curr;
        zprev = z0;
        
        wsol = w0;
        zsol = z0;
        while  wsol(end) >= wMin && wsol(end) <= wMax
            wpred = wsol(end) + pred.step*direction;
            if wpred > wMax || wpred < wMin
                break;
            end
            
            iPredict = max(length(wsol)-6,1):length(wsol);
            if length(iPredict) > 1
                zpred = interp1(wsol(iPredict),zsol(:,iPredict)',wpred,'pchip','extrap')';
            else
                zpred = zsol(:,end);
                zpred(end) = wpred;
            end
           
            %now try to solve
            xpred = zpred(1:end-1);
            Xpred(:,problem.iNL) = unpackdof(xpred,hbm.harm.NFreq-1,problem.NNL,hbm.harm.iRetainNL);
            sol = hbm_solve(hbm,problem,wpred,A,Xpred);
            sol.x = packdof(sol.X(:,problem.iNL),hbm.harm.iRetainNL);
            
            z = [sol.x; sol.w];
            t = z - zprev;
            corr.step = norm2(z - zprev);
            
            curr(end+1) = hbm_frf_results(z,t,pred,corr,hbm,problem);
            
            %unpack outputs           
            if ~any(isnan(curr(end).z))
                curr(end).flag = 'Success';
                pred.step = min(max(pred.step * hbm.cont.C,hbm.cont.min_step),hbm.cont.max_step);
                results(end+1) = curr(end);
                wsol(end+1) = results(end).w;
                zsol(:,end+1) = results(end).z;
                hbm_frf_plot('data',hbm,problem,results(end));
                prog.NFail = 0;
                if wsol(end) >= wMax || wsol(end) <= wMin
                    break;
                end
                zprev = curr(end).z;
            else
                curr(end).flag = 'Failed: No convergence';
                pred.step = pred.step * hbm.cont.c;
                prog.NFail = prog.NFail + 1;
                hbm_frf_plot('err',hbm,problem,curr(end));
                if prog.NFail > hbm.cont.maxfail
                    prog.Status = 'Too many failed iterations';
                    break;
                end
            end
        end
        
        %add on final point
        t = zEnd - zprev;
       
        pred.step = norm(zEnd - zprev);
        corr.step = norm(zEnd - zprev);
        corr.it = 0;
        curr(end+1) = hbm_frf_results(zEnd,t,pred,corr,hbm,problem);
        results(end+1) = curr(end);
        
        hbm_frf_plot('close',hbm,problem,[]);
        
    case 'predcorr'
       prog.NStep = 1;
       prog.NFail = 0;
       prog.NIter = 0;
                
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
        
        Zprev = z0./problem.Zscale;
        F = hbm_frf_constraints(Zprev,hbm,problem);
        J = hbm_frf_jacobian(Zprev,hbm,problem);
        Tprev = get_tangent(J);
        Tprev = Tprev * sign(Tprev(end)) * sign(wEnd - w0);
        
        Zend = zEnd./problem.Zscale;
        J = hbm_frf_jacobian(Zend,hbm,problem);
        Tend = get_tangent(J);
        Tend = Tend * sign(Tend(end)) * sign(wEnd - w0);
        
        pred.step = hbm.cont.step0;
        corr.step = pred.step;
        corr.it = 0;

        curr = hbm_frf_results(Zprev,Tprev,pred,corr,hbm,problem);
        results = curr;
        
        hbm_frf_plot('init',hbm,problem,init);
        fprintf('STEP    PRED    CORR   STATUS  INFO   ITER   TOT    FREQ      ')
        fprintf('Z(%d)       ',1:length(z0))
        fprintf('\n')
        fprintf('%3d   %6.4f  %6.4f    %s       %3s   %3d   %3d   %6.2f   ',prog.NStep,pred.step,corr.step,'S','Ini',corr.it,prog.NIter,results(end).w)
        fprintf('%+5.2e  ',results(end).z)
        fprintf('\n')       
        
        zsol = results.z;
        tsol = results.t;

        while norm(Zprev - Zend) > hbm.cont.max_step
            
            %predictor
            switch hbm.cont.predcorr.predictor
                case 'linear'
                    Zpred = Zprev + pred.step*Tprev;
                case 'quadratic'
                    if length(results) < 5
                        Zpred = Zprev + pred.step*Tprev;
                    else
                    	Zpred = polynomial_predictor(Zsol(:,end-2:end),Tsol(:,end-2:end),pred.step);
                    end
                case 'cubic'
                    if length(results) < 5
                        Zpred = Zprev + pred.step*Tprev;
                    else
                    	Zpred = polynomial_predictor(Zsol(:,end-3:end),Tsol(:,end-3:end),pred.step);
                    end
            end
            corr.it = 0;
            Zlast = Zpred + Inf;
            
            %corrector
            switch hbm.cont.predcorr.corrector
                case 'pseudo'
                    bConverged = 0;
                    Z = Zpred;
                    T = Tprev;

                    while corr.it <= hbm.cont.predcorr.maxit
                        Zlast = Z;
                        J = hbm_frf_jacobian(Z,hbm,problem);
                        F = hbm_frf_constraints(Z,hbm,problem);
                        if hbm.cont.predcorr.bMoorePenrose
                            Z = Zlast - J\F;
                        else
                            B = [J; T'];
                            R = [J*T; 0];
                            Q = [F; 0];
                            W = T - B\R;
                            T = normalise(W);
                            Z = Zlast - B\Q;
                        end
                        corr.it = corr.it + 1;
                        if ~(any(abs(Z - Zlast) > hbm.cont.xtol) || any(abs(F) > hbm.cont.ftol))
                            bConverged = 1;
                            break;
                        end
                    end
                    if hbm.cont.predcorr.bMoorePenrose
                        T = get_tangent(J);
                        T = sign(T'*Tprev)*T;
                    end
                case 'arclength'
                    corr.Zprev = Zprev;
                    corr.Tprev = Tprev;
                    corr.step = pred.step;
                    switch hbm.cont.predcorr.solver
                        case 'ipopt'
                            [Z,info] = fipopt('',Zpred,@hbm_arclength_constraints,ipopt_opt,hbm,problem,corr);
                            corr.it = info.iter;
                            bConverged = info.status == 0;
                        case 'fsolve'
                            [Z,F,status,out] = fsolve(@hbm_arclength_constraints,Zpred,fsolve_opt,hbm,problem,corr);
                            corr.it = out.iterations + 1;
                            bConverged = status == 1;
                    end
                    J = hbm_frf_jacobian(Z,hbm,problem);
                    T = get_tangent(J);
                    T = sign(T'*Tprev)*T;
            end
            corr.step = norm(Z - Zprev)*sign((Z-Zprev)'*Tprev);

            prog.NStep = prog.NStep + 1;
            prog.NIter = prog.NIter + corr.it;

            %prepare for plots
            curr(end+1) = hbm_frf_results(Z,T,pred,corr,hbm,problem);
            
            if bConverged && curr(end).sCorr >= hbm.cont.min_step && curr(end).sCorr <= hbm.cont.max_step && curr(end).w > 0
                %success
                status = 'S';
                hbm_frf_plot('data',hbm,problem,curr(end));
                
                if corr.it <= hbm.cont.num_iter_increase
                    pred.step = min(curr(end).sPred * hbm.cont.C,hbm.cont.max_step/curr(end).sCorr*curr(end).sPred);
                    curr(end).flag = 'Success: Increasing step size';
                    info = 'Inc';
                elseif corr.it >= hbm.cont.num_iter_reduce
                    pred.step = max(curr(end).sPred / hbm.cont.C,hbm.cont.min_step/curr(end).sCorr*curr(end).sPred);
                    curr(end).flag = 'Success: Reducing step size';
                    info = 'Red';
                else
                    curr(end).flag = 'Success';
                    info = '';
                end
                
                %store the data
                results(end+1) = curr(end);
                
                if bUpdateScaling
                    problem = hbm_scaling(problem,hbm,results(end));
                end
                
                zsol(:,end+1) = results(end).z;
                tsol(:,end+1) = results(end).t;
                Zsol = zsol./(repmat(problem.Zscale,1,size(zsol,2)));
                Tsol = normalise(tsol./(repmat(problem.Zscale,1,size(tsol,2))));

                Zend = zEnd./problem.Zscale;

                Zprev = Zsol(:,end);
                Tprev = Tsol(:,end);

                prog.NFail = 0;
            else
            	%failed
                hbm_frf_plot('err',hbm,problem,curr(end));
                status = 'F';
                
                if curr(end).w < 0
                    %gone to negative frequencies (somehow)
                    curr(end).flag = 'Failed: Negative frequency';
                    info = 'Neg';
                elseif corr.it <= hbm.cont.predcorr.maxit
                    %converge but step size is unacceptable
                    if curr(end).sCorr < 0
                        %backwards
                        pred.step = max(pred.step * hbm.cont.c,hbm.cont.min_step);
                        curr(end).flag = 'Failed: Wrong direction';
                        info = 'Bwd';
                    elseif curr(end).sCorr > hbm.cont.max_step
                        %too large
                        pred.step = max(0.9 * pred.step * hbm.cont.max_step / curr(end).sCorr,hbm.cont.min_step);
                        curr(end).flag = 'Failed: Step too large';
                        info = 'Lrg';
                    elseif curr(end).sCorr < hbm.cont.min_step
                        %too small
                        pred.step = min(1.1 * pred.step * hbm.cont.min_step / curr(end).sCorr,hbm.cont.max_step);
                        curr(end).flag = 'Failed: Step too small';
                        info = 'Sml';
                    else
                        pred.step = max(pred.step * hbm.cont.c,hbm.cont.min_step);
                        curr(end).flag = 'Failed: Other error';
                        info = 'Oth';
                    end
                else
                    %failed to converge
                    if any(abs(Z - Zlast) > hbm.cont.xtol)
                        curr(end).flag = 'Failed: No convergence';
                        info = 'xtl';
                    elseif any(abs(F) > hbm.cont.ftol)
                        curr(end).flag = 'Failed: Constraints violated';
                        info = 'ftl';
                    else
                        curr(end).flag = 'Failed';
                        info = '';
                    end
                    pred.step = max(pred.step * hbm.cont.c,hbm.cont.min_step);
                end
                
                prog.NFail = prog.NFail + 1;
                if prog.NFail > hbm.cont.maxfail
                    prog.Status = 'Too many failed iterations';
                    break
                elseif curr(end).w < 0
                    prog.Status = 'Negative frequency';
                    break
                end
            end
            
            fprintf('%3d   %6.4f  %6.4f    %s       %3s   %3d   %3d   %6.2f   ',prog.NStep,curr(end).sPred,curr(end).sCorr,status,info,corr.it,prog.NIter,curr(end).w)
            fprintf('%+5.2e  ',curr(end).z)
            fprintf('\n')
        end
        
        %add on final point
        pred.step = norm(Zend - Zprev);
        corr.step = norm(Zend - Zprev);
        corr.it = 0;
        curr(end+1) = hbm_frf_results(Zend,Tend,pred,corr,hbm,problem);
        results(end+1) = curr(end);

        hbm_frf_plot('close',hbm,problem,[]);
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
        prob = coco_add_func(prob, 'alg', @hbm_coco_constraints, @hbm_coco_jacobian, data, 'zero','u0', z0./problem.Zscale);
        prob = coco_add_pars(prob, 'pars', length(z0), 'w');
        prob = coco_add_slot(prob, 'plot', @hbm_coco_callback, data, 'bddat');
        
        prob = coco_set(prob,'cont','h0',hbm.cont.step0,'h_min',hbm.cont.min_step,'h_max',hbm.cont.max_step,'ItMX',hbm.cont.coco.ItMX,'NPR',hbm.cont.coco.NPR);
        prob = coco_set(prob,'ep','bifus',false);
        bd = coco(prob, name, [], 1, 'w', [wMin wMax]./problem.wscale);
        hbm_frf_plot('close',hbm,problem,[]);
        
        %extract the solutions
        lab_col = coco_bd_col(bd, 'TYPE');
        NPts = sum(cellfun(@(x)~isempty(x),lab_col));
        pred.step = 0;
        corr.step = 0;
        corr.it = 0;
        for i = 1:NPts
            chart = coco_read_solution(name,i,'chart');
            pred.step = chart.R;
            corr.step = chart.R;
            Z = chart.x(1:end-1);
            T = chart.t(1:end-1);
            results(i) = hbm_frf_results(Z,T,pred,corr,hbm,problem);
        end
        fclose all;
        try
            rmdir(['data' filesep name],'s')
        end
        cd(currdir);
        prog.Status = 'success';
end

NPts = length(results);

if ~isfield(results,'W')
    for i = 1:NPts
        results(i).W = (hbm.harm.kHarm*(hbm.harm.rFreqRatio.*hbm.harm.rFreqBase)')*results(i).w;
    end
end

if 0%~isfield(results,'L')
    for i = 1:NPts
        w0 = results(i).w*hbm.harm.rFreqRatio;
        results(i).L = hbm_floquet(hbm,problem,w0,results(i).U,results(i).X);
    end
end


%% Predictor
function X_extrap = polynomial_predictor(Z,dZ,s_extrap)
s = norm2(diff(Z,[],2));
s = cumsum([0 s]);
N = length(s);
if ~isempty(dZ)
    [A,Ad] = poly_mat(s,N);
    p = [A;Ad]\[Z dZ]';
else
    A = poly_mat(s,N);
    p = A\Z';
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

%% Constraints and Jacobian
function c = hbm_frf_constraints(Z,hbm,problem)
%unpack the inputs
x = Z(1:end-1).*problem.xscale;
w = Z(end).*problem.wscale;
A = problem.A;

w0 = w * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

c = hbm_balance3d('func',hbm,problem,w0,u,x);
c = c ./ problem.Fscale;

function J = hbm_frf_jacobian(Z,hbm,problem)
%unpack the inputs
x = Z(1:end-1).*problem.xscale;
w = Z(end).*problem.wscale;
A = problem.A;

w0 = w * hbm.harm.rFreqRatio;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

Jx = hbm_balance3d('jacob',hbm,problem,w0,u,x);
Dw = hbm_balance3d('derivW',hbm,problem,w0,u,x);

J = [Jx Dw];
J = J .* problem.Jscale;

%% Arclength files
function [c,J] = hbm_arclength_constraints(Z,hbm,problem,corr)
c = hbm_frf_constraints(Z,hbm,problem);
sgn = sign((Z - corr.Zprev)' * corr.Tprev);
s = norm(Z - corr.Zprev) * sgn;
c(end+1) = s - corr.step;
if nargout > 1
    J = hbm_arclength_jacobian(Z,hbm,problem,corr);
end

function J = hbm_arclength_jacobian(Z,hbm,problem,corr)
J = hbm_frf_jacobian(Z,hbm,problem);
sgn = sign((Z - corr.Zprev)' * corr.Tprev);
J(end+1,:) = sgn*((Z - corr.Zprev)'+eps)/(1*(norm(Z - corr.Zprev)+eps));

%% Coco functions
function [data, res] = hbm_coco_callback(prob, data, command, varargin)
hbm = data.hbm;
problem = data.problem;

switch command
    case 'init'
        x = prob.efunc.x0(1:end-1).*problem.xscale;
        init.X = unpackdof(x,hbm.harm.NHarm,problem.NDof,hbm.harm.iRetain);
        init.w = prob.efunc.x0(end).*problem.wscale;
        init.A = data.A;
        hbm_frf_plot('init',hbm,problem,init);
    case 'data'
        chart = varargin{1};
        x = chart.x(1:end-2).*problem.xscale;
        curr.X  = unpackdof(x,hbm.harm.NHarm,problem.NDof,hbm.harm.iRetain);
        curr.w = chart.x(end).*problem.wscale;
        curr.A = data.A;
        hbm_frf_plot('data',hbm,problem,curr);
end

res = {};

function [data,c] = hbm_coco_constraints(prob, data, u)
c = hbm_frf_constraints(u,data.hbm,data.problem);

function [data, J] = hbm_coco_jacobian(prob, data, u)
J = hbm_frf_jacobian(u,data.hbm,data.problem);

%% Utilities
function y = norm2(x)
y = sqrt(sum(x.^2,1));

function y = normalise(x)
scale = repmat(norm2(x),size(x,1),1);
y = x ./scale;

function t = get_tangent(J)
[U,S,V] = svd(J,0);
t = V(:,end);
t = t./norm(t);

function curr = hbm_frf_results(Z,tangent,pred,corr,hbm,problem)
w = Z(end).*problem.wscale;
xnl = Z(1:end-1).*problem.xscale;
Xnl = unpackdof(xnl,hbm.harm.NHarm,problem.NNL,hbm.harm.iRetainNL);
A = problem.A;
t = normalise(tangent.*problem.Zscale);

curr.z = [xnl; w];
curr.t = t;

curr.sCorr = corr.step;
curr.sPred = pred.step;
curr.it = corr.it;
curr.flag = '';

curr.w = w;
curr.U = A*feval(problem.excite,hbm,problem,curr.w*hbm.harm.rFreqRatio);
curr.X = hbm_recover(hbm,problem,curr.w*hbm.harm.rFreqRatio(1),curr.U,Xnl);
curr.F = hbm_output3d(hbm,problem,curr.w*hbm.harm.rFreqRatio,curr.U,curr.X);
curr.A = A;