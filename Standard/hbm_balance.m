function varargout = hbm_balance(command,hbm,problem,w,u,x)
NDofTot = hbm.harm.NComp*problem.NDof;
NNLTot = hbm.harm.NComp*problem.NNL;

ii = find(hbm.harm.NHarm ~= 0);
r = hbm.harm.rFreqRatio(ii);
w0 = w * r + hbm.harm.wFreq0(ii);

A  = (hbm.lin.Ak + w0*hbm.lin.Ac{ii} + w0^2*hbm.lin.Am{ii});
B  = (hbm.lin.Bk + w0*hbm.lin.Bc{ii} + w0^2*hbm.lin.Bm{ii});
dA0dw = r*(hbm.lin.Ac{ii} + 2*w0*hbm.lin.Am{ii});
dBdw  = r*(hbm.lin.Bc{ii} + 2*w0*hbm.lin.Bm{ii});
[A,R,dAdw,dRdw] = hbm_reduce(hbm,problem,A,dA0dw);

switch command
    case 'func' %F, used by hbm_frf & hbm_bb
        cl = B*u - hbm.lin.b - A*x;
        if hbm.bIncludeNL
            cnl = hbm_nonlinear('func',hbm,problem,w0,x,u);
        else
            cnl = 0*cl;
        end
        c = R*(cl - cnl);
        varargout{1} = c;
    case 'jacob' %dF_dX, used by hbm_frf & hbm_bb
        Jl = -A;
        if hbm.bIncludeNL
            if hbm.options.bAnalyticalDerivs
                [Jx,Jxdot,Jxddot] = hbm_nonlinear({'jacobX','jacobXdot','jacobXddot'},hbm,problem,w0,x,u);
                Jnl1 = Jx + w0*Jxdot + w0^2*Jxddot;
            else
                c0 = hbm_nonlinear('func',hbm,problem,w0,x,u);
                h = 1E-12;
                Jnl2 = zeros(NDofTot,NNLTot);
                x0 = x;
                u0 = u;
                for i = 1:NNLTot
                    x = x0;
                    x(i) = x(i) + h;
                    c = hbm_nonlinear('func',hbm,problem,w0,x,u0);
                    Jnl2(:,i) = (c-c0)./h;
                end
            end
            if hbm.options.bAnalyticalDerivs
                Jnl = Jnl1;
            else
                Jnl = Jnl2;
            end
        else
            Jnl = 0*Jl;
        end
        J = R*(Jl - Jnl);
        varargout{1} = J;
    case 'derivW' %dF_dw, used by hbm_frf & hbm_bb
        Dl = dBdw*u - dAdw*x;
        cl = B*u - hbm.lin.b - A*x;
        if hbm.bIncludeNL
            if hbm.dependence.xdot || hbm.dependence.w
                if hbm.options.bAnalyticalDerivs
                    [Jxdot,Jxddot,Judot,Juddot,Dw] = hbm_nonlinear({'jacobXdot','jacobXddot','jacobUdot','jacobUddot','derivW'},hbm,problem,w0,x,u);
                    Dxdot = r*Jxdot*x;
                    Dxddot = 2*r*w0*Jxddot*x;
                    Dudot = r*Judot*u;
                    Duddot = 2*r*w0*Juddot*u;
                    Dnl1 = Dxdot + Dxddot + Dudot + Duddot + Dw;
                else
                    cnl = hbm_nonlinear('func',hbm,problem,w0,x,u);
                    h = 1E-10;
                    c = hbm_nonlinear('func',hbm,problem,w0+h,x,u);
                    Dnl2 = r*(c-cnl)./h;
                end
                if hbm.options.bAnalyticalDerivs
                    Dnl = Dnl1;
                else
                    Dnl = Dnl2;
                end
            else
                Dnl = 0*Dl;
            end
        else
            Dnl = 0*Dl;
        end
        D = R*(Dl - Dnl) + dRdw*(cl-cnl);
        varargout{1} = D;
    case 'derivA' %dF_dA, used in hbm_bb
        Dl = B*u;
        if hbm.bIncludeNL
            [Ju,Judot,Juddot] = hbm_nonlinear({'jacobU','jacobUdot','jacobUddot'},hbm,problem,w0,x,u);
            Dnl = (Ju + w0*Judot + w0^2*Juddot)*u;
        else
            Dnl = 0*Dl;
        end
        D = R*(Dl - Dnl);
        varargout{1} = D;
    case 'floquet0'
        D0 = -hbm_balance('jacob',hbm,problem,w0,u,x);
        varargout{1} = D0;
    case 'floquet1'
        D1l = hbm.lin.floquet.D1xdot + 2*hbm.lin.floquet.D1xddot{ii}*w0;
        if hbm.bIncludeNL
            [D1xdot,D1xddot] = hbm_nonlinear({'floquet1xdot','floquet1xddot'},hbm,problem,w0,x,u);
            D1nl = D1xdot + 2*D1xddot*w0;
        else
            D1nl = 0*D1l;
        end
        varargout{1} = R*(D1l - D1nl);
    case 'floquet2'
        D2l = hbm.lin.floquet.D2;
        if hbm.bIncludeNL
            D2nl = hbm_nonlinear({'floquet2'},hbm,problem,w0,x,u);
        else
            D2nl = 0*D2l;
        end
        varargout{1} = R*(D2l - D2nl);
end