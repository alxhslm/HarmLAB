function results = hbm_bb(hbm,problem,A0,w0,X0,AEnd,wEnd,XEnd)
NDof = problem.NDof;

%first solve @ A0
sol = hbm_resonance(hbm,problem,w0,A0,X0);
X0 = sol.X;
w0 = sol.w0;
if any(isnan(abs(X0(:))))
    error('Failed to solve initial problem')
end
u0 = packdof(A0*feval(problem.excite,hbm,problem,w0*hbm.harm.rFreqRatio));
x0 = packdof(X0);
H0 = hbm_objective('func',hbm,problem,w0*hbm.harm.rFreqRatio,x0,u0);

sol = hbm_resonance(hbm,problem,wEnd,AEnd,XEnd);
XEnd = sol.X;
wEnd = sol.w0;
if any(isnan(abs(XEnd(:))))
    error('Failed to solve final problem')
end
xEnd = packdof(XEnd);
uEnd = packdof(AEnd*feval(problem.excite,hbm,problem,wEnd*hbm.harm.rFreqRatio));
HEnd = hbm_objective('func',hbm,problem,wEnd*hbm.harm.rFreqRatio,xEnd,uEnd);

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
    problem = hbm_scaling(problem,hbm,x0,w0,A0);
    bUpdateScaling = 1;
end
 
x0 = packdof(X0);
xEnd = packdof(XEnd);

AMax = max(AEnd,A0);
AMin = min(AEnd,A0);

problem.AMax = AMax;
problem.AMin = AMin;
problem.A0   = A0;
problem.AEnd = AEnd;

switch hbm.cont.method
      case 'eitherside'
        [fig,ax] = hbm_bb_plot('init',hbm,problem,A0,H0,w0);
        A = [A0 AEnd];
        Ascale = mean([A0 AEnd]);
        x = [x0 xEnd];
        w = [w0 wEnd];
        u = [u0 uEnd];       
        H = [H0 HEnd];
        
        hbm_bb_plot('data',hbm,problem,A0,H0,w0);
        hbm_bb_plot('data',hbm,problem,AEnd,HEnd,wEnd);
        step = hbm.cont.step0;
        Aupper  = AEnd;
        Alower  = A0;
        iLower = 1;
        iUpper = 1;
        num_err = 0;
        bUpper = 0;
        
        order = 3;
        Aplot = linspace(A0,AEnd,1000);
        B =  poly_mat(log(Aplot),order+1);
        for i = 1:2
            hPred(i) = plot(ax(i),Aplot,Aplot+NaN,'b-');
        end
        while true
            if bUpper
                Apred = Aupper - step*Ascale;
                if Apred < Alower
                    break
                end
            else
                Apred = Alower + step*Ascale;
                if Apred > Aupper
                    break
                end
            end

            if iUpper > 1 && iLower >  1
                xpred = polynomial_predictor(log(A),x,log(Apred),order);
                [wpred,pw] = polynomial_predictor(log(A),w,log(Apred),order);
                [Hpred,ph] = polynomial_predictor(log(A),H,log(Apred),order);
                set(hPred(1),'YData',(B*pw)');
                set(hPred(2),'YData',(B*ph)');
                drawnow
            elseif bUpper
                xpred = xEnd;
                wpred = wEnd;
            else
                xpred = x0;
                wpred = w0;
            end
            
            Xpred = unpackdof(xpred,hbm.harm.NFreq-1,problem.NDof);
            [Xsol,wsol] = hbm_resonance(hbm,problem,wpred,Apred,Xpred);
            
            xsol = packdof(Xsol);
            usol = packdof(Apred*feval(problem.excite,hbm,problem,wsol*hbm.harm.rFreqRatio));
            Hsol = hbm_objective('func',hbm,problem,wsol*hbm.harm.rFreqRatio,xsol,usol);
            
            if ~any(isnan(xsol))
                step = step * hbm.cont.C;
                
                x = [x(:,1:iLower) xsol x(:,end-iUpper+1:end)];
                u = [u(:,1:iLower) usol u(:,end-iUpper+1:end)];
                w = [w(1:iLower) wsol w(end-iUpper+1:end)];
                H = [H(1:iLower) Hsol H(end-iUpper+1:end)];             
                A = [A(1:iLower) Apred A(end-iUpper+1:end)];  
                
                if bUpper
                    Aupper = Apred;
                    iUpper = iUpper + 1;
                else
                    Alower = Apred;
                    iLower = iLower + 1;
                end
                
                bUpper = ~bUpper;
                
                hbm_bb_plot('data',hbm,problem,Apred,Hsol,wsol);
                num_err = 0;
            else
                step = step * hbm.cont.c;
                bUpper = ~bUpper;
                
                num_err = num_err + 1;
                hbm_bb_plot('err',hbm,problem,Apred,Hsol,wsol);
                if num_err > hbm.cont.maxfail
                    break;
                end
            end
        end
        hbm_bb_plot('close',hbm,problem,[],0,0);
    case 'none'
        [~,ax] = hbm_bb_plot('init',hbm,problem,A0,H0,w0);
        A = A0;
        Ascale = mean([A0 AEnd]);
        x = x0;
        w = w0;
        u = packdof(A*feval(problem.excite,hbm,problem,w*hbm.harm.rFreqRatio));
        h = hbm_objective('func',hbm,problem,w*hbm.harm.rFreqRatio,x,u);
        
        hbm_bb_plot('data',hbm,problem,A0,H0,w0);
        step = hbm.cont.step0;
        direction = sign(AEnd-A0)*Ascale;

        order = 3;
        Aplot = logspace(0,0.07,100).^sign(AEnd-A0);
        B =  poly_mat(log(Aplot),order+1);
        for i = 1:2
            hPred(i) = plot(ax(i),Aplot*A0,Aplot+NaN,'b--');
        end
        
        while A(end) <= AMax && A(end) >= AMin
            Apred = A(end) + step*direction;
            if Apred > AMax || Apred < AMin
                break;
            end
            
            iPredict = max(length(A)-6,1):length(A);
            Aoff = A(iPredict(1));
            xpred = polynomial_predictor(A(iPredict),x(:,iPredict),Apred);
            [wpred,pw] = polynomial_predictor(A(iPredict),w(iPredict),Apred);
            [Hpred,ph] = polynomial_predictor(A(iPredict),h(iPredict),Apred);
            
            set(hPred(1),'XData',Aplot*Aoff,'YData',(B(:,1:length(pw))*pw)');
            set(hPred(2),'XData',Aplot*Aoff,'YData',(B(:,1:length(ph))*ph)');
            drawnow
            
            Xpred = unpackdof(xpred,hbm.harm.NFreq-1,problem.NDof);
            [Xsol,wsol] = hbm_resonance(hbm,problem,wpred,Apred,Xpred);
            xsol = packdof(Xsol);
            usol = packdof(Apred*feval(problem.excite,hbm,problem,wsol*hbm.harm.rFreqRatio));
            Hsol = hbm_objective('func',hbm,problem,wsol*hbm.harm.rFreqRatio,xsol,usol);
            if ~any(isnan(xsol))
                step = min(max(step * hbm.cont.C,hbm.cont.min_step),hbm.cont.max_step);
                x(:,end+1) = xsol;
                u(:,end+1) = usol;
                w(end+1) = wsol;
                h(end+1) = Hsol;               
                A(end+1) = Apred;
                hbm_bb_plot('data',hbm,problem,Apred,Hsol,wsol);
                num_err = 0;
                if A(end) >= AMax || A(end) <= AMin
                    break;
                end
            else
                step = step * hbm.cont.c;
                num_err = num_err + 1;
                hbm_bb_plot('err',hbm,problem,Apred,Hsol,wsol);
                if num_err > hbm.cont.maxfail
                    break;
                end
            end
        end
        if num_err > hbm.cont.maxfail
            x(:,end+1) = NaN;
            w(end+1) = NaN;
            A(end+1) = NaN;
        else
            x(:,end+1) = xEnd;
            w(end+1)   = wEnd;
            A(end+1)   = AEnd;
        end
        hbm_bb_plot('close',hbm,problem,[],0,0);
    case 'predcorr'
        w = w0;
        x = x0;
        A = A0;
        h = H0;
        Xprev = [x0; w0; A0]./problem.Xscale;
        F0 = hbm_pseudo_constraints(Xprev,hbm,problem);
        J0 = hbm_pseudo_jacobian(Xprev,hbm,problem);
                
        hbm_bb_plot('init',hbm,problem,A,h,w);
        fprintf('STEP    PRED    CORR   STATUS  INFO     FREQ        AMP     ITER   TOT      ')
        fprintf('X(%d)       ',1:length(x0))
        fprintf('\n')       
        
        fprintf('%3d   %6.4f  %6.4f    %s      %3s    %6.2f      %6.2f    %3d   %3d   ',0,0,0,'S','   ',w0,A0,0,0)
        fprintf('%+5.2e  ',x0)
        fprintf('\n')
        
        step = hbm.cont.step0;
        num_fail = 0;
        num_step = 0;
        num_iter_tot = 0;
        flag = {};
        
        tangent_prev = pseudo_null(J0);
        tangent_prev = tangent_prev * sign(tangent_prev(end)) * sign(AEnd - A0);
        t0 = tangent_prev.*problem.Xscale;
        t = t0;
        
        while A(end) <= AMax && A(end) >= AMin
            
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
                            F = hbm_arclength_constraints(X,hbm,problem,A,Xprev,tangent_prev,step);
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
            aCurr = X(end).*problem.Xscale(end);
            wCurr = X(end-1).*problem.Xscale(end-1);
            xCurr = X(1:end-2).*problem.Xscale(1:end-2);
            tCurr = tangent.*problem.Xscale;
            uCurr = packdof(aCurr*feval(problem.excite,hbm,problem,wCurr*hbm.harm.rFreqRatio));
            HCurr = hbm_objective('func',hbm,problem,wCurr*hbm.harm.rFreqRatio,xCurr,uCurr);
            
            step_prev = norm(X - Xprev);

            if bConverged && sign((X - Xprev)'*tangent_prev) > 0
                %success
                status = 'S';
                
                if num_iter <= hbm.cont.num_iter_increase
                    step = min(step * hbm.cont.C,hbm.cont.max_step/step_prev*step);
                    flag{end+1} = 'Success: Increasing step size';
                    info = 'Inc';
                elseif num_iter >= hbm.cont.num_iter_reduce
                    step = max(step / hbm.cont.C,hbm.cont.min_step/step_prev*step);
                    flag{end+1} = 'Success: Reducing step size';
                    info = 'Red';
                else
                    flag{end+1} = 'Success';
                    info = '';
                end
                
                %just right
                num_fail = 0;
                Xprev = X;

                if num_iter < hbm.cont.num_iter_reduce
                    step = min(step * hbm.cont.C,step * hbm.cont.max_step / step_prev);
                elseif num_iter > hbm.cont.num_iter_increase
                    step = step / hbm.cont.C;
                end

                tangent_prev = sign(tangent'*tangent_prev)*tangent;

                %store the data
                w(end+1) = wCurr;
                A(end+1) = aCurr;
                x(:,end+1) = xCurr;
                h(end+1) = HCurr;
                t(:,end+1) = tCurr;

                if bUpdateScaling
                    problem = hbm_scaling(problem,hbm,x(:,end),w(end),A(end));
                end

                Xsol = [x;w;A]./(repmat(problem.Xscale,1,size(x,2))); 
                tsol = normalise(t./(repmat(problem.Xscale,1,size(t,2))));

                hbm_bb_plot('data',hbm,problem,aCurr,HCurr,wCurr);

            else
                %failed
                status = 'F';
                
                if wCurr(end) < 0
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
                
                hbm_bb_plot('err',hbm,problem,aCurr,HCurr,wCurr);
                num_fail = num_fail + 1;
                if num_fail > 4
                    break
                end
            end
            fprintf('%3d   %6.4f  %6.4f    %s      %3s    %6.2f      %6.2f    %3d   %3d   ',num_step,step,step_prev,status,info,wCurr(end),aCurr(end),num_iter,num_iter_tot)
            fprintf('%+5.2e  ',xCurr(:,end))
            fprintf('\n')
        end
        hbm_bb_plot('close',hbm,problem,[],[],0);
    case 'coco'
        currdir = pwd;
        cd(desktoproot);
        data.hbm = hbm;
        data.problem = problem;
        
        prob = coco_prob();
        prob = coco_add_func(prob, 'alg', @hbm_coco_constraints, @hbm_coco_jacobian, data, 'zero','u0', [x0(:); w0; A0]./Xscale);
        prob = coco_add_pars(prob, 'pars', length(x0)+2, 'A');
        prob = coco_add_slot(prob, 'plot', @hbm_coco_callback, data, 'bddat');

        prob = coco_set(prob,'cont','h0',hbm.cont.step0,'h_min',hbm.cont.min_step,'h_max',hbm.cont.max_step,'ItMX',hbm.cont.coco.ItMX,'NPR',hbm.cont.coco.NPR);
        bd = coco(prob, 'temp', [], 1, 'A', [AMin AMax]./Xscale(end));
        hbm_bb_plot('close',hbm,problem,[],[],[]);
        
        %extract the solutions
        lab_col = coco_bd_col(bd, 'TYPE');
        NPts = sum(cellfun(@(x)~isempty(x),lab_col));
        x = zeros(hbm.harm.NComp*NDof,NPts);
        w = zeros(NPts,1);
        A = zeros(NPts,1);
        for i = 1:NPts
            chart = coco_read_solution('temp',i,'chart');
            x(:,i) = chart.x(1:end-3);
            w(i) = chart.x(end-2);
            A(i) = chart.x(end);
        end
        w = w .* Xscale(end-1);
        A = A .* Xscale(end);
        x = x .* (Xscale(1:end-2)*(0*w'+1));
        fclose all;
        try
            rmdir('data','s')
        end
        cd(currdir);
end

X = unpackdof(x,hbm.harm.NFreq-1,NDof);

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


if ~exist('h','var')
    h = 0*w;
    for i = 1:NPts
        h(i) = hbm_objective('complex',hbm,problem,w(i)*hbm.harm.rFreqRatio,x(:,i),u(:,i));
    end
end

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
    'H',h,...
    'L',lambda);

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

% function [X_extrap,p] = polynomial_predictor(s,X,s_extrap,order)
% N = length(s);
% if nargin < 4
%     order = N-1;
% end
% A = poly_mat(s,order+1);
% p = A\X';
% B = poly_mat(s_extrap,order+1);
% X_extrap = (B*p)';
% 
% function [A,Ad] = poly_mat(s,N)
% A = zeros(length(s),N);
% Ad = zeros(length(s),N);
% for i = 1:N
%     A(:,i) = s.^(i-1);
%     if nargout > 1
%         if i > 1
%             Ad(:,i) = (i-1)*s.^(i-2);
%         else
%             Ad(:,i) = 0*s;
%         end
%     end
% end

function t = pseudo_null(A)
[U,S,V] = svd(A,0);
t = V(:,end);
t = t./norm(t);

function c = hbm_pseudo_constraints(X,hbm,problem)
%unpack the inputs
Xscale = problem.Xscale;
x = X(1:end-2).*Xscale(1:end-2);
w = X(end-1).*Xscale(end-1);
w0 = w * hbm.harm.rFreqRatio;
A = X(end).*Xscale(end);

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

c = hbm_balance3d('func',hbm,problem,w0,u,x);
c(end+1) = resonance_condition(hbm,problem,w,x,A);

function J = hbm_pseudo_jacobian(X,hbm,problem)

%unpack the inputs
Xscale = problem.Xscale;
x = X(1:end-2).*Xscale(1:end-2);
w = X(end-1).*Xscale(end-1);
w0 = w * hbm.harm.rFreqRatio;
A = X(end).*Xscale(end);

% X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof);
U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

Jx = hbm_balance3d('jacob',hbm,problem,w0,u,x);
Dw = hbm_balance3d('derivW',hbm,problem,w0,u,x);
Da = hbm_balance3d('derivA',hbm,problem,w0,u,x);

res0 = resonance_condition(hbm,problem,w,x,A);
h = 1E-10;
x0 = x;
for i = 1:length(x)
    x = x0;
    x(i) = x(i) + h;
    drdx(i) = (resonance_condition(hbm,problem,w,x,A) - res0)/h;
    i
end

drdw = (resonance_condition(hbm,problem,w+h,x0,A) - res0)/h;
drdA = (resonance_condition(hbm,problem,w,x0,A+h) - res0)/h;

J = [Jx   Dw    Da;
     drdx drdw drdA];

J = J .* repmat(Xscale(:)',size(J,1),1);

function Hw = resonance_condition(hbm,problem,w,x,A)
w0 = w * hbm.harm.rFreqRatio;
U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);
Jx = hbm_balance3d('jacob',hbm,problem,w0,u,x);
Dw = hbm_balance3d('derivW',hbm,problem,w0,u,x);

Hw = hbm_peak(x,w0,A,hbm,problem,Jx,Dw);

function dHdw = hbm_peak(x,w0,A,hbm,problem,Jx,Dw)
x0 = x;
obj0 = hbm_obj(x0,w0,A,hbm,problem);

h = 1E-10;
Gx = 0*x;
for i = 1:length(x)
    x = x0;
    x(i) = x(i) + h;
    Gx(i) = (hbm_obj(x,w0,A,hbm,problem)-obj0)/h;
end

Gw = (hbm_obj(x0,w0+h,A,hbm,problem)-obj0)/h;

dHdw = Gw - Gx'*(Jx\Dw);

function obj = hbm_obj(x,w0,A,hbm,problem)
X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof);
U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

f = hbm_output3d(hbm,problem,w0,u,x);
F = unpackdof(f,hbm.harm.NFreq-1,problem.NOutput);

obj = feval(problem.obj,X,U,F,hbm,problem,w0);

%% Coco functions

function [data res] = hbm_coco_callback(prob, data, command, varargin)
hbm = data.hbm;
problem = data.problem;
Xscale = problem.Xscale;

switch command
    case 'init'
        res = {};
        hbm_bb_plot('init',hbm,problem,[],0,0);
    case 'data'
        chart = varargin{1};
        w = chart.x(end-2).*Xscale(end-1);
        a = chart.x(end).*Xscale(end);
        hbm_bb_plot('data',hbm,problem,a,w);
        res = {};
end

function [data,c] = hbm_coco_constraints(prob, data, u)
hbm = data.hbm;
problem = data.problem;
Xscale = problem.Xscale;

x = u(1:end-2).*Xscale(1:end-2);
w = u(end-1).*Xscale(end-1);
w0 = w * hbm.harm.rFreqRatio;
A = u(end).*Xscale(end);

X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof);
U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

c = hbm_balance3d('func',hbm,problem,w0,u,x);

c(end+1) = resonance_condition(hbm,problem,w,x,A);

function [data, J] = hbm_coco_jacobian(prob, data, u)
hbm = data.hbm;
problem = data.problem;
Xscale = problem.Xscale;

x = u(1:end-2).*Xscale(1:end-2);
w = u(end-1).*Xscale(end-1);
w0 = w * hbm.harm.rFreqRatio;
A = u(end).*Xscale(end);

X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof);
U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

Jx = hbm_balance3d('jacob',hbm,problem,w0,u,x);
Dw = hbm_balance3d('derivW',hbm,problem,w0,u,x);
Da = hbm_balance3d('derivA',hbm,problem,w0,u,x)/A;

%now get the derivatives of the resonance condition
[~,Hwx, Hw2, HwA] = resonance_condition(hbm,problem,w,x,A);

J = [Jx  Dw  Da;
     Hwx Hw2 HwA];

J = J .* repmat(Xscale(:)',size(J,1),1);

%% Arclength files

function f = hbm_arclength_constraints(X,hbm,problem,Xprev,tangent_prev,step)
Xscale = problem.Xscale;
x = X(1:end-2).*Xscale(1:end-2);
w = X(end-1).*Xscale(end-1);
w0 = w * hbm.harm.rFreqRatio;
A = X(end).*Xscale(end);
Xcurr = X;

X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof);
U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

c = hbm_balance3d('func',hbm,problem,w0,u,x);

sgn = sign((X - Xprev)' * tangent_prev);
s = norm(X - Xprev) * sgn;

c(end+1) = resonance_condition(hbm,problem,w,x,A);

f = [c;
     s - step];
      
function J = hbm_arclength_jacobian(X,hbm,problem,Xprev,tangent_prev,step)
Xscale = problem.Xscale;
Xcurr = X;

x = X(1:end-2).*Xscale(1:end-2);
w = X(end-1).*Xscale(end-1);
w0 = w * hbm.harm.rFreqRatio;
A = X(end).*Xscale(end);

X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof);
U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

sgn = sign((X - Xprev)' * tangent_prev);


Jx = hbm_balance3d('jacob',hbm,problem,w0,u,x);
Dw = hbm_balance3d('derivW',hbm,problem,w0,u,x);
Da = hbm_balance3d('derivA',hbm,problem,w0,u,x)/A;

%now get the derivatives of the resonance condition
[~,Hwx, Hw2, HwA] = resonance_condition(hbm,problem,w,x,A);

J = [Jx  Dw  Da;
     Hwx Hw2 HwA];

J = J .* repmat(Xscale(:)',size(J,1),1);

J(end+1,:) = sgn*((Xcurr - Xprev)'+eps)./(norm(Xcurr - Xprev)+eps);
J = sparse(J);