function varargout = hbm_objective(part,hbm,problem,w,xnl,u)

[x,R] = hbm_recover(hbm,problem,w,u,xnl);

NOutput = problem.res.NOutput;
NInput = problem.res.NInput;

r = hbm.harm.rFreqRatio;
w0 = w .* r + hbm.harm.wFreq0;

if problem.res.iHarm > 1
    iRe = (1:NOutput)' + 2*(problem.res.iHarm - 2)*NOutput + NOutput;
    iIm = (1:NOutput)' + 2*(problem.res.iHarm - 2)*NOutput + 2*NOutput;
    
    jRe = (1:NInput)' + 2*(problem.res.iHarm - 2)*NInput + NInput;
    jIm = (1:NInput)' + 2*(problem.res.iHarm - 2)*NInput + 2*NInput;
else
    iRe = (1:NOutput)';
    iIm = [];
    
    jRe = (1:NInput)';
    jIm = [];
end

kk = hbm.harm.kHarm*(hbm.harm.rFreqBase.*hbm.harm.rFreqRatio)';
D = [0 -1;
     1  0];
D = blkdiag(0,kron(diag(kk(2:end)),D));

%% output
Wx = kron(D,eye(problem.NDof));
switch problem.res.output
    case 'none'
        Fnl = 0*x;
    case 'fnl'
        Fnl = hbm_nonlinear3d({'func'},hbm,problem,w0,xnl,u);
    case 'x'
        Fnl = x;
    case 'xdot'
        Fnl = w*Wx*x;
    case 'xddot'
        Fnl = w^2*Wx*Wx*x;
end

Fb = Fnl(iRe);
if ~isempty(iIm)
    Fb = Fb + 1i*Fnl(iIm);
end
Fb = problem.res.ROutput*Fb;

%% input
Wu = kron(D,eye(problem.NInput));
switch problem.res.input
    case 'unity'
        Fex = 0*u+1;
    case 'fe'
        Ju = prod(w0)*hbm.lin.Bx + hbm.lin.Bk;
        for k = 1:2
            Ju = Ju + (w0(k)*hbm.lin.Bc{k} + w0(k)^2*hbm.lin.Bm{k});
        end
        Fex = Ju*u;
    case 'u'
        Fex = u;
    case 'udot'
        Fex = w * Wu*u;
    case 'uddot'
        Fex = w^2 * Wu*Wu*u;
end

Fe = Fex(jRe);
if ~isempty(jIm)
    Fe = Fe + 1i*Fex(jIm);
end
Fe = problem.res.RInput*Fe;

if ~iscell(part)
    part = {part};
end
varargout = cell(1,length(part));

for i = 1:length(part)
    switch part{i}
        case 'complex'
            H = Fb./Fe;
            varargout{i} = H;
        case 'func'
            H = abs(Fb./Fe);
            varargout{i} = H;
        case 'jacobX'
            
            %nl
            switch problem.res.output
                case 'none'
                    Jx = 0*Wx;
                case 'fnl'
                    Jx = hbm_nonlinear3d({'jacobX'},hbm,problem,w0,xnl,u);
                case 'x'
                    Jx = eye(problem.NDof*hbm.harm.NComp);
                case 'xdot'
                    Jx = Wx;
                case 'xddot'
                    Jx = Wx*Wx;
            end
            
            if ~isempty(iIm)
                dFbdx = (real(Fb)*Jx(iRe,:) + imag(Fb)*Jx(iIm,:))./(abs(Fb) + eps);
            else
                dFbdx = (real(Fb)*Jx(iRe,:))./(abs(Fb) + eps);
            end
            dFbdx = problem.res.ROutput*dFbdx;

            %excitation
            dFedx = 0;
            
            %put it all together
            dHdx = (dFbdx.*abs(Fe) - abs(Fb).*dFedx)./abs(Fe).^2;
            
            varargout{i} = dHdx*R;
        case 'derivW'
           
            %nl
            switch problem.res.output
                case 'none'
                    Dw_nl = 0*x;
                case 'fnl'
                    Dw = hbm_nonlinear3d({'derivW'},hbm,problem,w0,xnl,u);
                    if ~iscell(Dw)
                        Dw = {Dw,0*Dw};
                    end
                    Dw_nl = 0*Dw{1};
                    for k = 1:2
                        Dw_nl = Dw_nl + r(k)*Dw{k};
                    end
                case 'x'
                    Dw_nl = zeros(problem.NDof*hbm.harm.NComp,1);
                case 'xdot'
                    Dw_nl = (Wx*x);
                case 'xddot'
                    Dw_nl = 2*w*(Wx*Wx*x);
            end
            if ~isempty(iIm)
                dFbdw = (real(Fb)*Dw_nl(iRe) + imag(Fb)*Dw_nl(iIm))./(abs(Fb) + eps);
            else
                dFbdw = (real(Fb)*Dw_nl(iRe))./(abs(Fb) + eps);
            end
            dFbdw = problem.res.ROutput*dFbdw;
            
            %excitation
            switch problem.res.input
                case 'unity'
                    Dw_u = 0*u;
                case 'fe'
                    Dw_u = (r(1)*w0(2) + r(2)*w0(1))*hbm.lin.Bx*u;
                    for k = 1:2
                        Dw_u = Dw_u + r(k)*(hbm.lin.Bc{k} + 2*w0(k)*hbm.lin.Bm{k})*u;
                    end
                case 'u'
                    Dw_u = zeros(problem.NInput*hbm.harm.NComp,1);
                case 'udot'
                    Dw_u = (Wu*u);
                case 'uddot'
                    Dw_u = 2*w*(Wu*Wu*u);
            end
            if ~isempty(iIm)
                dFedw = (real(Fe)*Dw_u(jRe) + imag(Fe)*Dw_u(jIm))./abs(Fe);
            else
                dFedw = (real(Fe)*Dw_u(jRe))./abs(Fe);
            end
            dFedw = problem.res.RInput*dFedw;
            
            dHdw = (dFbdw.*abs(Fe) - abs(Fb).*dFedw)./abs(Fe).^2;
            
            varargout{i} = dHdw;
    end
end