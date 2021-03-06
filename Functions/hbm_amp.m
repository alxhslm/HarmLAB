function [results,curr] = hbm_amp(hbm,problem,w0,A0,X0,AEnd,XEnd)
problem.type = 'amp';
problem.w0 = w0;

%first solve @ A0
sol = hbm_solve(hbm,problem,w0,A0,X0);
x0 = packdof(sol.X);
if any(isnan(abs(x0(:))))
    error('Failed to solve initial problem')
end
z0 = [x0;A0];

init.X = sol.X;
init.w = w0;
init.A = A0;

sol = hbm_solve(hbm,problem,w0,AEnd,XEnd);
xEnd = packdof(sol.X);
if any(isnan(abs(xEnd(:))))
    error('Failed to solve final problem')
end
zEnd = [xEnd;AEnd];

hbm.bIncludeNL = 1;

if isfield(problem,'xscale')
    xscale = [problem.xscale'; repmat(problem.xscale',hbm.harm.NFreq-1,1)*(1+1i)];   
    problem.Xscale = packdof(xscale)*sqrt(length(xscale));
    problem.Ascale = mean([A0 AEnd]);
    problem.Fscale = problem.Xscale*0+1;
    problem.Zscale = [problem.Xscale; problem.Ascale];
    problem.Jscale = (1./problem.Fscale(:))*problem.Zscale(:)';
    bUpdateScaling = 0;
else
    problem = hbm_scaling(problem,hbm,init);
    bUpdateScaling = 1;
end
 
AMax = max(AEnd,A0);
AMin = min(AEnd,A0);

problem.AMax = AMax;
problem.AMin = AMin;
problem.A0   = A0;
problem.AEnd = AEnd;


prog.Status = 'success';

switch hbm.cont.method
    case 'none'
        hbm_amp_plot('init',hbm,problem,init);
        Ascale = mean([A0 AEnd]);
        
        pred.step = hbm.cont.step0;
        direction = sign(AEnd-A0)*Ascale;
        
        problem.Xscale = 0*problem.Xscale + 1;
        problem.Ascale = 1;
        problem.Zscale = 0*problem.Zscale + 1;
        
        z = z0;
        J = hbm_amp_jacobian(z,hbm,problem);
        t = get_tangent(J);
        
        corr.step = 0;
        corr.it = 0;
        
        curr = hbm_amp_results(z,t,pred,corr,hbm,problem);
        results = curr;
        zprev = z0;
        
        Asol = A0;
        zsol = z0;
        
        while Asol(end) <= AMax && Asol(end) >= AMin
            Apred = Asol(end) + pred.step*direction;
            if Apred > AMax || Apred < AMin
                break;
            end
            
            iPredict = max(length(Asol)-6,1):length(Asol);
            if length(iPredict) > 1
                zpred = interp1(Asol(iPredict),zsol(:,iPredict)',Apred,'pchip','extrap')';
            else
                zpred = zsol(:,end);
                zpred(end) = Apred;
            end
            
            %now try to solve
            xpred = zpred(1:end-1);
            Xpred = unpackdof(xpred,hbm.harm.NFreq-1,problem.NDof);
            sol = hbm_solve(hbm,problem,w0,Apred,Xpred);
            sol.x = packdof(sol.X);
            
            z = [sol.x; sol.A];
            t = z - zprev;
            corr.step = norm2(z - zprev);
            
            curr(end+1) = hbm_amp_results(z,t,pred,corr,hbm,problem);
            
            if ~any(isnan(curr(end).z))
                curr(end).flag = 'Success';
                pred.step = min(max(pred.step * hbm.cont.C,hbm.cont.min_step),hbm.cont.max_step);
                results(end+1) = curr(end);
                Asol(end+1) = results(end).A;
                zsol(:,end+1) = results(end).z;
                hbm_amp_plot('data',hbm,problem,results(end));
                prog.NFail = 0;
                if Asol(end) >= AMax || Asol(end) <= AMin
                    break;
                end
                zprev = curr(end).z;
            else
                curr(end).flag = 'Fail';
                pred.step = pred.step * hbm.cont.c;
                prog.NFail = prog.NFail + 1;
                hbm_amp_plot('err',hbm,problem,curr(end));
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
        curr(end+1) = hbm_amp_results(zEnd,t,pred,corr,hbm,problem);
        results(end+1) = curr(end);
        
        hbm_amp_plot('close',hbm,problem,[]);
        
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
        F = hbm_amp_constraints(Zprev,hbm,problem);
        J = hbm_amp_jacobian(Zprev,hbm,problem);
        Tprev = get_tangent(J);
        Tprev = Tprev * sign(Tprev(end)) * sign(AEnd - A0);

        Zend = zEnd./problem.Zscale;
        J = hbm_amp_jacobian(Zend,hbm,problem);
        Tend = get_tangent(J);
        Tend = Tend * sign(Tend(end)) * sign(AEnd - A0);
        
        pred.step = hbm.cont.step0;
        corr.step = pred.step;
        corr.it = 0;

        curr = hbm_amp_results(Zprev,Tprev,pred,corr,hbm,problem);
        results = curr;
        
        hbm_amp_plot('init',hbm,problem,init);
        fprintf('STEP    PRED    CORR   STATUS  INFO   ITER   TOT     AMP     ')
        fprintf('Z(%d)       ',1:length(z0))
        fprintf('\n')       
        
        fprintf('%3d   %6.4f  %6.4f    %s      %3s   %3d    %3d    %6.2f   ',prog.NStep,pred.step,corr.step,'S','Ini',corr.it,prog.NIter,results(end).A)
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
                        J = hbm_amp_jacobian(Z,hbm,problem);
                        F = hbm_amp_constraints(Z,hbm,problem);
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
                    J = hbm_amp_jacobian(Z,hbm,problem);
                    T = get_tangent(J);
                    T = sign(T'*Tprev)*T;
            end
            corr.step = norm(Z - Zprev)*sign((Z-Zprev)'*Tprev);

            prog.NStep = prog.NStep + 1;
            prog.NIter = prog.NIter + corr.it;

            %prepare for plots
            curr(end+1) = hbm_amp_results(Z,T,pred,corr,hbm,problem);
            
            if  bConverged && curr(end).sCorr >= hbm.cont.min_step && curr(end).sCorr <= hbm.cont.max_step && curr(end).A > 0
                %success
                status = 'S';
                hbm_amp_plot('data',hbm,problem,curr(end));
                
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
                hbm_amp_plot('err',hbm,problem,curr(end));
                status = 'F';
                
                if curr(end).A < 0
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
                elseif curr(end).A < 0
                    prog.Status = 'Negative frequency';
                    break
                end
            end
            fprintf('%3d   %6.4f  %6.4f    %s      %3s   %3d    %3d    %6.2f  ',prog.NStep,curr(end).sPred,curr(end).sCorr,status,info,corr.it,prog.NIter,curr(end).A)
            fprintf('%+5.2e  ',curr(end).z)
            fprintf('\n')
        end
        
        %add on final point
        pred.step = norm(Zend - Zprev);
        corr.step = norm(Zend - Zprev);
        corr.it = 0;
        curr(end+1) = hbm_amp_results(Zend,Tend,pred,corr,hbm,problem);
        results(end+1) = curr(end);

        hbm_amp_plot('close',hbm,problem,[]);
end

NPts = length(results);

if ~isfield(results,'W')
    for i = 1:NPts
        w0 = results(i).w*hbm.harm.rFreqRatio + hbm.harm.wFreq0;
        results(i).W = hbm.harm.kHarm*(hbm.harm.rFreqBase.*w0)';
    end
end

if ~isfield(results,'L')
    results = hbm_floquet(hbm,problem,results);
end

results = hbm_excitation_forces(problem,results);

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
function c = hbm_amp_constraints(Z,hbm,problem)
%unpack the inputs
x = Z(1:end-1).*problem.Xscale;
A = Z(end).*problem.Ascale;
w = problem.w0;

w0 = w * hbm.harm.rFreqRatio + hbm.harm.wFreq0;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

c = hbm_balance3d('func',hbm,problem,w,u,x);
c = c ./ problem.Fscale;

function J = hbm_amp_jacobian(Z,hbm,problem)
%unpack the inputs
x = Z(1:end-1).*problem.Xscale;
A = Z(end).*problem.Ascale;
w = problem.w0;

w0 = w * hbm.harm.rFreqRatio + hbm.harm.wFreq0;

U = A*feval(problem.excite,hbm,problem,w0);
u = packdof(U);

Jx = hbm_balance3d('jacob' ,hbm,problem,w,u,x);
Da = hbm_balance3d('derivA',hbm,problem,w,u,x);

J = [Jx  Da];
J = J .* problem.Jscale;

%% Arclength files
function [c,J] = hbm_arclength_constraints(Z,hbm,problem,corr)
c = hbm_amp_constraints(Z,hbm,problem);
sgn = sign((Z - corr.Zprev)' * corr.Tprev);
s = norm(Z - corr.Zprev) * sgn;
c(end+1) = s - corr.step;
if nargout > 1
    J = hbm_arclength_jacobian(Z,hbm,problem,corr);
end

function J = hbm_arclength_jacobian(Z,hbm,problem,corr)
J = hbm_amp_jacobian(Z,hbm,problem);
sgn = sign((Z - corr.Zprev)' * corr.Tprev);
J(end+1,:) = sgn*((Z - corr.Zprev)'+eps)/(1*(norm(Z - corr.Zprev)+eps));

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

function curr = hbm_amp_results(Z,tangent,pred,corr,hbm,problem)
A = Z(end).*problem.Ascale;
x = Z(1:end-1).*problem.Xscale;
w = problem.w0;
t = normalise(tangent.*problem.Zscale);

curr.z = [x; A];
curr.t = t;

curr.sCorr = corr.step;
curr.sPred = pred.step;
curr.it = corr.it;
curr.flag = '';

w0 = w*hbm.harm.rFreqRatio + hbm.harm.wFreq0;

curr.w = w;
curr.X = unpackdof(x,hbm.harm.NHarm,problem.NDof);
curr.U = A*feval(problem.excite,hbm,problem,w0);
curr.F = hbm_output3d(hbm,problem,curr.w,curr.U,curr.X);
curr.A = A;