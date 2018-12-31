function results = hbm_frf(hbm,problem,A,w0,x0,wEnd,xEnd)
NDof = problem.NDof;

%first solve @ w0
sol = hbm_solve(hbm,problem,w0,A,x0);
if any(isnan(abs(x0(:))))
    error('Failed to solve initial problem')
end
x0 = packdof(sol.X,hbm.harm.iRetain);
u0 = packdof(sol.U);
f0 = packdof(sol.F);

sol = hbm_solve(hbm,problem,wEnd,A,xEnd);
xEnd = packdof(sol.X,hbm.harm.iRetain);

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
    problem = update_scaling(problem,hbm,x0,w0);
    bUpdateScaling = 1;
end

wMax = max(wEnd,w0);
wMin = min(wEnd,w0);

problem.wMax = wMax;
problem.wMin = wMin;
problem.w0   = w0;
problem.wEnd = wEnd;

err = 'success';

switch hbm.options.cont_method
    case 'none'
        w = w0;
        x = x0;
        s = [];
        it = [];
        hbm_frf_plot('init',hbm,problem,hbm.cont.step0,x,w,A);
        step = hbm.cont.step0;
        while  w(end) >= wMin && w(end) <= wMax
            if ~any(isnan(x(:,end)))
                x0 = x(:,end);
            end
            wCurr = w(end) + step*problem.wscale*sign(wEnd-w0);
            X = hbm_solve(hbm,problem,wCurr,A,unpackdof(x0,hbm.harm.NFreq-1,problem.NDof,hbm.harm.iRetain));
            if any(isnan(X))
                err = {'failed to solve'};
                break;
            end
            if bUpdateScaling
                problem.x0scale = abs(X(1,:));
                problem.xhscale = max(abs(X(2:end,:)),[],1);
            end
            x(:,end+1) = packdof(X);
            w(end+1) = wCurr;
            s(end+1) = step;
            it(end+1) = 1;
            hbm_frf_plot('data',hbm,problem,x(:,end),w(end),A);
        end
        hbm_frf_plot('close',hbm,problem,[],[],[],[]);
        debug = struct();
    case 'pseudo'
        fun_constr   = @(X)hbm_pseudo_constraints(X,hbm,problem,A);
        fun_jacobian = @(X)hbm_pseudo_jacobian(X,hbm,problem,A);
        
        w = w0; wCurr = w0;
        x = x0; xCurr = x0;
        u = u0; uCurr = u0;
        f = f0; fCurr = f0;
        s = []; sCorrCurr = []; sPredCurr = [];
        it = []; itCurr = [];
        flag = {};
%         lambda = lambda0;
        
        Xprev = [x0; w0]./problem.Xscale;

        num_step = 1;
        num_iter = 0;
        num_fail = 0;
        num_iter_tot = 0;
        
        step = hbm.cont.step0;
        
        hbm_frf_plot('init',hbm,problem,x,w,A);
        fprintf('STEP    PRED    CORR   STATUS  INFO  ITER   TOT   FREQ       ')
        fprintf('X(%d)       ',1:length(x0))
        fprintf('\n')
        fprintf('%3d   %6.4f  %6.4f    %s       %3s  %3d   %3d   %6.2f   ',num_step,step,step,'S','   ',num_iter,num_iter_tot,w)
        fprintf('%+5.2e  ',x)
        fprintf('\n')       
        
        F = hbm_pseudo_constraints(Xprev,hbm,problem,A);
        J = hbm_pseudo_jacobian(Xprev,hbm,problem,A);
        tangent_prev = pseudo_null(J);
        tangent_prev = tangent_prev * sign(tangent_prev(end)) * sign(wEnd - w0);
        t0 = tangent_prev.*problem.Xscale;
        XCurr = Xprev;
        t = t0; tCurr = t0;
        Xend = [xEnd;wEnd]./problem.Xscale;
        
        while norm(Xprev - Xend) > hbm.cont.max_step
            
            if length(w) < 5
                %linear
                X = Xprev + step*tangent_prev;
            else
                %cubic
                X = polynomial_predictor(Xsol(:,end-3:end),tsol(:,end-3:end),step);
            end
            num_iter = 0;
            Xlast = X + Inf;
            tangent = tangent_prev;
            
            while (any(abs(X - Xlast) > hbm.cont.xtol) || any(abs(F) > hbm.cont.ftol)) && num_iter <= hbm.cont.pseudo.maxit
                Xlast = X;
                J = hbm_pseudo_jacobian(X,hbm,problem,A);
                F = hbm_pseudo_constraints(X,hbm,problem,A);
                if hbm.cont.pseudo.bMoorePenrose
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
            end
%             X(end) = abs(X(end));
            
            num_step = num_step + 1;
            num_iter_tot = num_iter_tot + num_iter;
            
            if hbm.cont.pseudo.bMoorePenrose
                tangent = pseudo_null(J);
                tangent = sign(tangent'*tangent_prev)*tangent;
            end
            
            if tangent(end) > 0
%                  1
            end
            
            %prepare for plots
            XCurr(:,end+1) = X;
            wCurr(end+1) = X(end).*problem.wscale;
            xCurr(:,end+1) = X(1:end-1).*problem.xscale;
            tCurr(:,end+1) = tangent.*problem.Xscale;
            uCurr(:,end+1) = packdof(A*feval(problem.excite,hbm,problem,wCurr(end)*hbm.harm.rFreqRatio));
            fCurr(:,end+1) = hbm_output3d(hbm,problem,wCurr(end)*hbm.harm.rFreqRatio,uCurr(:,end),xCurr(:,end));
            sCorrCurr(end+1) = norm(X - Xprev)*sign((X-Xprev)'*tangent_prev);
            sPredCurr(end+1) = step;
            itCurr(end+1) = num_iter;
            
            if num_iter <= hbm.cont.pseudo.maxit && sCorrCurr(end) >= hbm.cont.min_step && sCorrCurr(end) <= hbm.cont.max_step && wCurr(end) > 0
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
                f(:,end+1) = fCurr(:,end);
                t(:,end+1) = tCurr(:,end);
                s(end+1)   = sCorrCurr(end);
                it(end+1)  = itCurr(end);
%                 lambda(:,end+1) = floquetMultipliers(hbm,problem,wCurr(end)*hbm.harm.rFreqRatio,uCurr(:,end),xCurr(:,end));
                
                if bUpdateScaling
                    problem = update_scaling(problem,hbm,x(:,end),w(end));
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
                elseif num_iter <= hbm.cont.pseudo.maxit
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
        XCurr(:,end+1) = Xend;
        wCurr(end+1)   = wEnd;
        xCurr(:,end+1) = xEnd;
        tCurr(:,end+1) = tCurr(:,end);
        uCurr(:,end+1) = packdof(A*feval(problem.excite,hbm,problem,wCurr(end)*hbm.harm.rFreqRatio));
        fCurr(:,end+1) = hbm_output3d(hbm,problem,wCurr(end)*hbm.harm.rFreqRatio,uCurr(:,end),xCurr(:,end));
        sCorrCurr(end+1) = norm(Xend - Xprev);
        sPredCurr(end+1) = norm(Xend - Xprev);
        flag{end+1} = 'Success';
        itCurr(end+1) = 0;
        
%         %floquet multipliers
%         lambda(:,end+1) = floquetMultipliers(hbm,problem,wCurr(end)*hbm.harm.rFreqRatio,uCurr(:,end),xCurr(:,end));

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
    case 'arclength'
        w = w0; wCurr = w0;
        x = x0; xCurr = x0;
        u = u0; uCurr = u0;
        f = f0; fCurr = f0;
        s = []; sCorrCurr = []; sPredCurr = [];
        it = []; itCurr = [];
        flag = {};
%         lambda = lambda0;
        
        problem = update_scaling(problem,hbm,x0,w0);
        Xprev = [x0; w0]./problem.Xscale;
        Xend = [xEnd;wEnd]./problem.Xscale;
        step = hbm.cont.step0;
        num_iter_tot = 0;
        
        XCurr = Xprev;
        
        hbm_frf_plot('init',hbm,problem,x,w,A);
        fprintf('STEP LAB    PAR(1)         ')
        fprintf('U(%d)          ',1:length(x0))
        fprintf('\n')
        fprintf('%3d %3d  %+0.5e   ',0,1,w)
        fprintf('%+0.5e   ',x)
        fprintf('\n')
        
        Jstr = [hbm.sparsity 0*x0+1;
                0*x0'+1 1];
        
        switch hbm.cont.arclength.solver
            case 'ipopt'
                ipopt_opt.print_level = 0;
                ipopt_opt.maxit = hbm.cont.arclength.maxit;
                ipopt_opt.ftol = hbm.cont.ftol;
                ipopt_opt.xtol = hbm.cont.xtol;
                ipopt_opt.jacob = @hbm_arclength_jacobian;
                ipopt_opt.jacobstructure = Jstr;
            case 'fsolve'
                fsolve_opt = optimoptions('fsolve',...
                                            'Display','off',...
                                            'TolFun',hbm.cont.ftol,...
                                            'TolX',hbm.cont.xtol,...
                                            'Jacobian','on',...
                                            'JacobPattern',Jstr,...
                                            'MaxIter',hbm.cont.arclength.maxit);
        end
        
        tangent_prev = 0*Xprev; tangent_prev(end) = sign(wEnd-w0);
        F = hbm_arclength_constraints(Xprev,hbm,problem,A,Xprev,tangent_prev,step);
        J = hbm_arclength_jacobian(Xprev,hbm,problem,A,Xprev,tangent_prev,step);
        tangent_prev = pseudo_null(J(1:end-1,:));
        tangent_prev = tangent_prev * sign(tangent_prev(end)) * sign(wEnd - w0);
        
        while norm(Xprev - Xend) > hbm.cont.max_step
            bConverged = 0;
            num_fail = 0;
            
            while ~bConverged               
                if length(wCurr) < 6
                    %linear
                    X0 = Xprev + step*tangent_prev;
                else
                    %cubic
                    X0 = polynomial_predictor(Xsol(:,end-3:end),[],step);
                end
                
                switch hbm.cont.arclength.solver
                    case 'ipopt'
                        [X,info] = fipopt([],X0,@hbm_arclength_constraints,ipopt,hbm,problem,A,Xprev,tangent_prev,step);
                        iter = info.iter;
                        bConverged = info.status == 0;
                        F = hbm_arclength_constraints(X,hbm,problem,A,Xprev,tangent_prev,step);
                    case 'fsolve'
                        [X,F,status,out] = fsolve(@hbm_arclength_constraints,X0,fsolve_opt,hbm,problem,A,Xprev,tangent_prev,step);
                        iter = out.iterations + 1;
                        bConverged = status == 1;
                end
                num_iter_tot = num_iter_tot + iter;
                
                %prepare for plots
                XCurr(:,end+1) = X;
                wCurr(end+1) = X(end).*problem.wscale;
                xCurr(:,end+1) = X(1:end-1).*problem.xscale;
                uCurr(:,end+1) = packdof(A*feval(problem.excite,hbm,problem,wCurr(end)*hbm.harm.rFreqRatio));
                fCurr(:,end+1) = hbm_output3d(hbm,problem,wCurr(end)*hbm.harm.rFreqRatio,uCurr(:,end),xCurr(:,end));
                sCorrCurr(end+1) = norm(X - Xprev);
                sPredCurr(end+1) = step;
                itCurr(end+1) = iter;
                
                if bConverged
                    if sCorrCurr(end) > hbm.cont.max_step
                        %too large step
                        step = 0.9 * step * hbm.cont.max_step / sCorrCurr(end);
                        flag{end+1} = 'Step too large';
                    elseif sCorrCurr(end) < hbm.cont.min_step
                        step = 1.1 * step * hbm.cont.min_step / sCorrCurr(end);
                        flag{end+1} = 'Step too small';
                    else
                        %just right
                        if iter <= hbm.cont.num_iter_increase
                            step = min(step * hbm.cont.C, step * hbm.cont.max_step / sCorrCurr(end));
                            flag{end+1} = 'Increasing step size';
                        elseif iter >= hbm.cont.num_iter_reduce
                            step = step / hbm.cont.C;
                            flag{end+1} = 'Reducing step size';
                        else
                            flag{end+1} = 'Success';
                        end

                        J = hbm_arclength_jacobian(X,hbm,problem,A,Xprev,tangent_prev,step);
                        hz = null(J(1:end-1,:));
                        tangent = sign(hz'*tangent_prev)*hz;
                        tangent = tangent .* problem.Xscale;
                        
%                         %floquet multipliers
%                         lambda(:,end+1) = floquetMultipliers(hbm,problem,wCurr(end)*hbm.harm.rFreqRatio,uCurr(:,end),xCurr(:,end));
                        
                        %store the data
                        w(end+1) = wCurr(end);
                        x(:,end+1) = xCurr(:,end);
                        u(:,end+1) = uCurr(:,end);
                        f(:,end+1) = fCurr(:,end);
                        s(end+1) = sCorrCurr(end);
                        
                        problem = update_scaling(problem,hbm,x(:,end),w(end));

                        Xprev = [x(:,end); w(end)]./ problem.Xscale;
                        tangent_prev = tangent ./ problem.Xscale;
                        Xsol = [x;w]./(repmat(problem.Xscale,1,size(x,2))); 
                        
                        hbm_frf_plot('data',hbm,problem,xCurr(:,end),wCurr(end),A);
                        
                        fprintf('%3d %3d  %+0.5e   ',num_iter_tot,length(w),w(end))
                        fprintf('%+0.5e   ',x(:,end))
                        fprintf('\n')
                    end
                else
                    %failed
                    flag{end+1} = 'Failed';
                    step = step*hbm.cont.c;
                    hbm_frf_plot('err',hbm,problem,xCurr(:,end),wCurr(end),A);
                    num_fail = num_fail + 1;
                    if num_fail > 4
                        break;
                    end
                end
                drawnow
            end
        end

        %add on final point
        XCurr(:,end+1) = Xend;
        wCurr(end+1)   = wEnd;
        xCurr(:,end+1) = xEnd;
        uCurr(:,end+1) = packdof(A*feval(problem.excite,hbm,problem,wCurr(end)*hbm.harm.rFreqRatio));
        fCurr(:,end+1) = hbm_output3d(hbm,problem,wCurr(end)*hbm.harm.rFreqRatio,uCurr(:,end),xCurr(:,end));
        sCorrCurr(end+1) = norm(Xend - Xprev);
        sPredCurr(end+1) = norm(Xend - Xprev);
        flag{end+1} = 'Success';
        itCurr(end+1) = 0;

        hbm_frf_plot('close',hbm,problem,[],[],[],[]);
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
        hbm_frf_plot('close',hbm,problem,[],[],[],[]);
        
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
        x = x .* (xscale*(0*w+1));
        fclose all;
        try
            rmdir(['data' filesep name],'s')
        end
        cd(currdir);
        debug = struct();
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

function problem = update_scaling(problem,hbm,x,w)
X = unpackdof(x,hbm.harm.NFreq-1,problem.NDof,hbm.harm.iRetain);
xdc  = max(abs(X(1,:)),1E-6);
xmax = max(max(abs(X(2:end,:)),[],1),1E-6);
xharm = repmat(xmax,hbm.harm.NFreq-1,1);
% xharm = max(abs(X(2:end,:)),1E-10);
xscale = [xdc; xharm*(1+1i)];
problem.xscale = packdof(xscale,hbm.harm.iRetain);
problem.wscale = w;
problem.Fscale = x*0+1;
problem.Xscale = [problem.xscale; problem.wscale];
problem.Jscale = (1./problem.Fscale(:))*problem.Xscale(:)';

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
    Jx = hbm_balance3d('jacob',hbm,problem,w0,u,x);
    Dw = hbm_balance3d('derivW',hbm,problem,w0,u,x);
    
    J = [Jx Dw];
    J = J .*problem.Jscale;
    J(end+1,:) = sgn*((X - Xprev)'+eps)/(1*(norm(X - Xprev)+eps)); %length(X)
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