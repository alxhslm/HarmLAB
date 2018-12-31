function varargout = hbm_derivatives(fun,var,t,x,xdot,xddot,u,udot,uddot,xalg,f,hbm,problem,w0)
NPts = size(x,1);

if ~iscell(var)
    var = {var};
end

for i = 1:length(var)
    switch var{i}
        case 'alg'
            df_dxalg = feval(problem.model,[fun '_alg'],t',x.',xdot.',xddot.',u.',udot.',uddot.',xalg.',hbm,problem,w0);
            if isempty(df_dxalg)
                if hbm.dependence.xalg
                    f0 = f.';
                    h = 1E-10;
                    df_dxalg = zeros(size(f,2),problem.NAlg,NPts);
                    for j = 1:problem.NAlg
                        xAlg2 = xalg.';
                        xAlg2(j,:) = xAlg2(j,:) + h;
                        f2 = feval(problem.model,fun,t',x.',xdot.',xddot.',u.',udot.',uddot.',xAlg2,hbm,problem,w0);
                        df_dxalg(:,j,:) = permute((f2-f0)./h,[1 3 2]);
                    end
                else
                    df_dxalg = zeros(size(f,2),problem.NAlg,NPts);
                end
            end
            varargout{i} = df_dxalg;
        case 'x'
            df_dx = feval(problem.model,[fun '_x'],t',x.',xdot.',xddot.',u.',udot.',uddot.',xalg.',hbm,problem,w0);
            if isempty(df_dx)
                if hbm.dependence.x
                    f0 = f.';
                    h = 1E-10;
                    df_dx = zeros(size(f,2),problem.NDof,NPts);
                    for j = 1:problem.NDof
                        x2 = x.';
                        x2(j,:) = x2(j,:) + h;
                        f2 = feval(problem.model,fun,t',x2,xdot.',xddot.',u.',udot.',uddot.',xalg.',hbm,problem,w0);
                        df_dx(:,j,:) = permute((f2-f0)./h,[1 3 2]);
                    end
                else
                    df_dx = zeros(size(f,2),problem.NDof,NPts);
                end
            end
            varargout{i} = df_dx;
        case 'xdot'
            df_dxdot = feval(problem.model,[fun '_xdot'],t',x.',xdot.',xddot.',u.',udot.',uddot.',xalg.',hbm,problem,w0);
            if isempty(df_dxdot)
                if hbm.dependence.xdot
                    f0 = f.';
                    h = 1E-10;
                    df_dxdot = zeros(size(f,2),problem.NDof,NPts);
                    for j = 1:problem.NDof
                        xdot2 =  xdot.';
                        xdot2(j,:) = xdot2(j,:) + h;
                        f2 = feval(problem.model,fun,t',x.',xdot2,xddot.',u.',udot.',uddot.',xalg.',hbm,problem,w0);
                        df_dxdot(:,j,:) = permute((f2-f0)./h,[1 3 2]);
                    end
                else
                    df_dxdot = zeros(size(f,2),problem.NDof,NPts);
                end
            end
            varargout{i} = df_dxdot;
        case 'xddot'
            df_dxddot = feval(problem.model,[fun '_xddot'],t',x.',xdot.',xddot.',u.',udot.',uddot.',xalg.',hbm,problem,w0);
            if isempty(df_dxddot)
                if hbm.dependence.xddot
                    f0 = f.';
                    h = 1E-10;
                    df_dxddot = zeros(size(f,2),problem.NDof,NPts);
                    for j = 1:problem.NDof
                        xddot2 =  xddot.';
                        xddot2(j,:) = xddot2(j,:) + h;
                        f2 = feval(problem.model,fun,t',x.',xdot.',xddot2,u.',udot.',uddot.',xalg.',hbm,problem,w0);
                        df_dxddot(:,j,:) = permute((f2-f0)./h,[1 3 2]);
                    end
                else
                    df_dxddot = zeros(size(f,2),problem.NDof,NPts);
                end
            end
            varargout{i} = df_dxddot;
        case 'u'
            df_du = feval(problem.model,[fun '_u'],t',x.',xdot.',xddot.',u.',udot.',uddot.',xalg.',hbm,problem,w0);
            if isempty(df_du)
                if hbm.dependence.u
                    f0 = f.';
                    h = 1E-10;
                    df_du = zeros(size(f,2),problem.NInput,NPts);
                    for j = 1:problem.NInput
                        u2 = u.';
                        u2(j,:) = u2(j,:) + h;
                        f2 = feval(problem.model,fun,t',x.',xdot.',xddot.',u2,udot.',uddot.',xalg.',hbm,problem,w0);
                        df_du(:,j,:) = permute((f2-f0)./h,[1 3 2]);
                    end
                else
                    df_du = zeros(size(f,2),problem.NInput,NPts);
                end
            end
            varargout{i} = df_du;
        case 'udot'
            df_dudot = feval(problem.model,[fun '_udot'],t',x.',xdot.',xddot.',u.',udot.',uddot.',xalg.',hbm,problem,w0);
            if isempty(df_dudot)
                if hbm.dependence.udot
                    f0 = f.';
                    h = 1E-10;
                    df_dudot = zeros(size(f,2),problem.NInput,NPts);
                    for j = 1:problem.NInput
                        udot2 =  udot.';
                        udot2(j,:) = udot2(j,:) + h;
                        f2 = feval(problem.model,fun,t',x.',xdot.',xddot.',u.',udot2,uddot.',xalg.',hbm,problem,w0);
                        df_dudot(:,j,:) = permute((f2-f0)./h,[1 3 2]);
                    end
                else
                    df_dudot = zeros(size(f,2),problem.NInput,NPts);
                end
            end
            varargout{i} = df_dudot;
         case 'uddot'
            df_duddot = feval(problem.model,[fun '_uddot'],t',x.',xdot.',xddot.',u.',udot.',uddot.',xalg.',hbm,problem,w0);
            if isempty(df_duddot)
                if hbm.dependence.udot
                    f0 = f.';
                    h = 1E-10;
                    df_duddot = zeros(size(f,2),problem.NInput,NPts);
                    for j = 1:problem.NInput
                        uddot2 =  uddot.';
                        uddot2(j,:) = uddot2(j,:) + h;
                        f2 = feval(problem.model,fun,t',x.',xdot.',xddot.',u.',udot.',uddot2,xalg.',hbm,problem,w0);
                        df_duddot(:,j,:) = permute((f2-f0)./h,[1 3 2]);
                    end
                else
                    df_duddot = zeros(size(f,2),problem.NInput,NPts);
                end
            end
            varargout{i} = df_duddot;
        case 'w'
            df_dw = feval(problem.model,[fun '_w'],t', x.',xdot.',xddot.',u.',udot.',uddot.',xalg.',hbm,problem,w0);
            if isempty(df_dw)
                f0 = f.';
                h = 1E-6;
                NFreq = length(w0);
                
                if hbm.dependence.w
                     for k = 1:length(w0)
                        w2 = w0;
                        t2 = t';
                        w2(k) = w2(k) + h;
                        t2(k,:) = w0(k)*t2(k,:)/w2(k);
                        f = feval(problem.model,fun,t2, x.',xdot.',xddot.',u.',udot.',uddot.',xalg.',hbm,problem,w2);
                        df_dw{k} = (f-f0)./h;
                    end
                else
                    df_dw = repmat({zeros(size(f,2),NPts)},1,NFreq);
                end
            end
            for k = 1:length(w0)
                df_dw{k} = df_dw{k}.';
            end
            varargout{i} = df_dw;
    end
end