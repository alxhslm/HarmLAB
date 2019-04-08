function varargout = hbm_balance(command,hbm,problem,w0,u,x)
NDofTot = hbm.harm.NComp*problem.NDof;

switch command
    case 'func' %F, used by hbm_frf & hbm_bb
        Jx = (hbm.lin.Ak{1} + w0*hbm.lin.Ac{1} + w0^2*hbm.lin.Am{1});
        Ju = (hbm.lin.Bk{1} + w0*hbm.lin.Bc{1} + w0^2*hbm.lin.Bm{1});
        cl = Ju*u - hbm.lin.b - Jx*x;
        if hbm.bIncludeNL
            cnl = hbm_nonlinear('func',hbm,problem,w0,x,u);
        else
            cnl = 0*cl;
        end
        c = cl - cnl;
        varargout{1} = c;
    case 'jacob' %dF_dX, used by hbm_frf & hbm_bb
        Jl = -(hbm.lin.Ak{1} + w0*hbm.lin.Ac{1} + w0^2*hbm.lin.Am{1});        
        if hbm.bIncludeNL
            if hbm.options.bAnalyticalDerivs
                [Jx,Jxdot,Jxddot] = hbm_nonlinear({'jacobX','jacobXdot','jacobXddot'},hbm,problem,w0,x,u);
                Jnl1 = Jx + w0*Jxdot + w0^2*Jxddot;
            else
                c0 = hbm_nonlinear('func',hbm,problem,w0,x,u);
                h = 1E-12;
                Jnl2 = zeros(NDofTot);
                x0 = x;
                u0 = u;
                for i = 1:NDofTot
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
        
        J = Jl - Jnl;
        varargout{1} = J;
    case 'derivW' %dF_dw, used by hbm_frf & hbm_bb
        Dl = (hbm.lin.Bc{1} + 2*w0*hbm.lin.Bm{1})*u - (hbm.lin.Ac{1} + 2*w0*hbm.lin.Am{1})*x;
        if hbm.bIncludeNL
            if hbm.options.bAnalyticalDerivs
                [Jxdot,Jxddot,Judot,Juddot,Dw] = hbm_nonlinear({'jacobXdot','jacobXddot','jacobUdot','jacobUddot','derivW'},hbm,problem,w0,x,u);
                Dxdot = Jxdot*x;
                Dxddot = 2*w0*Jxddot*x;
                Du = (Judot + 2*w0*Juddot)*u;
                Dnl1 = Dxdot + Dxddot + Du + Dw;
            else
                c0 = hbm_nonlinear('func',hbm,problem,w0,x,u);
                h = 1E-10;
                c = hbm_nonlinear('func',hbm,problem,w0+h,x,u);
                Dnl2 = (c-c0)./h;
            end
            if hbm.options.bAnalyticalDerivs
                Dnl = Dnl1;
            else
                Dnl = Dnl2;
            end
        else
            Dnl = 0*Dl;
        end
        D = Dl - Dnl;
        varargout{1} = D;
    case 'derivA' %dF_dA, used in hbm_bb
        Dl = (hbm.lin.Bk{1} + w0*hbm.lin.Bc{1} + w0^2*hbm.lin.Bm{1})*u;
        if hbm.bIncludeNL
            [Ju,Judot,Juddot] = hbm_nonlinear({'jacobU','jacobUdot','jacobUddot'},hbm,problem,w0,x,u);
            Dnl = (Ju + w0*Judot + w0^2*Juddot)*u;
        else
            Dnl = 0*Dl;
        end
        D = Dl - Dnl;
        varargout{1} = D;
    case 'floquet0'
        D0 = -hbm_balance('jacob',hbm,problem,w0,u,x);
        varargout{1} = D0;
    case 'floquet1'
        D1l = hbm.lin.floquet.D1xdot + 2*hbm.lin.floquet.D1xddot{1}*w0;
        if hbm.bIncludeNL
            [D1xdot,D1xddot] = hbm_nonlinear({'floquet1xdot','floquet1xddot'},hbm,problem,w0,x,u);
            D1nl = D1xdot + 2*D1xddot*w0;
        else
            D1nl = 0*D1l;
        end
        varargout{1} = D1l - D1nl;
    case 'floquet2'
        D2l = hbm.lin.floquet.D2;
        if hbm.bIncludeNL
            D2nl = hbm_nonlinear({'floquet2'},hbm,problem,w0,x,u);
        else
            D2nl = 0*D2l;
        end
        varargout{1} = D2l - D2nl;
        %% All cases below are currently unused - they are required for computing the resonance, but we now use FD
    case 'jacobW' %d2F_dwdX, to be used in hbm_bb
        Dl = - (hbm.lin.Ac{1} + 2*w0*hbm.lin.Am{1});
        if hbm.bIncludeNL
            Jxdot = hbm_nonlinear('jacobXdot',hbm,problem,w0,x,u);
            Dnl = Jxdot;
        else
            Dnl = 0*Dl;
        end
        D = Dnl + Dl;
        varargout{1} = D;
    case 'derivW2' %d2F_dw2, to be used in hbm_bb
        Dl = 2*hbm.lin.Bm{1}*u - 2*hbm.lin.Am{1}*x;
        Dnl = 0*Dl;
        D = Dnl + Dl;
        varargout{1} = D;
    case 'derivWA' %d2F_dwdA, to be used in hbm_bb
        Dl = (hbm.lin.Bc{1} + 2*w0*hbm.lin.Bm{1})*u;
        if hbm.bIncludeNL
            Judot = hbm_nonlinear({'jacobUdot'},hbm,problem,w0,x,u);
            Dnl = Judot*u;
        else
            Dnl = 0*Dl;
        end
        D = Dnl + Dl;
        varargout{1} = D;
end