function varargout = hbm_nonlinear3d(command,hbm,problem,w0,Xp,Up)
if hbm.options.bUseStandardHBM
    [varargout{1:nargout}] = hbm_nonlinear(command,hbm,problem,w0(1),Xp,Up);
    return;
end

NInput = problem.NInput;
NDof = problem.NDof;

NFreq = hbm.harm.NFreq;
NComp = hbm.harm.NComp;
NDofTot = NComp*NDof;
NInputTot = NComp*NInput;

%work out the time domain
X = unpackdof(Xp,NFreq-1,NDof);
U = unpackdof(Up,NFreq-1,NInput);

%get the time series
States = hbm_states3d(w0,X,U,hbm);

%push through the nl system
States.f = feval(problem.model,'nl',States,hbm,problem);

%finally convert into the frequency domain
F = hbm.nonlin.FFT*States.f.';
Fp = packdof(F);

if ~iscell(command)
    command = {command};
end

varargout = cell(1,length(command));

ijacobx = hbm.nonlin.hbm.ijacobx;
ijacobu = hbm.nonlin.hbm.ijacobu;

for o = 1:length(command)
    switch command{o}
        case 'func'
            %all done
            varargout{o} = Fp;
        case 'derivW' %df_dw
            if ~hbm.dependence.w
                Dw = repmat({zeros(NDofTot,1)},1,2);
            else
                dfhbm_dw = hbm_derivatives('nl','w',States,hbm,problem);
                if isfield(problem,'derivW')
                    Dw = feval(problem.derivW,dfhbm_dw,hbm,problem);
                else
                    for n = 1:2
                        Dw{n} = packdof(hbm.nonlin.FFT*dfhbm_dw{n});
                    end
                end
            end
            varargout{o} = Dw;
        case 'jacobX'  %df_dX = Jx*X
            if ~hbm.dependence.x
                Jx = zeros(NDofTot,NDofTot);
            else
                States.df_dx = hbm_derivatives('nl','x',States,hbm,problem);
                
                if isfield(problem,'jacobX')
                    Jx = feval(problem.jacobX,States,hbm,problem);
                else
                    Jx = sum(hbm.nonlin.hbm.Jx.*States.df_dx(ijacobx,ijacobx,:),3);
                end
            end
            varargout{o} = Jx;
        case 'jacobU' %df_dU = Ju*U
            if ~hbm.dependence.u
                Ju = zeros(NDofTot,NInputTot);
            else
                States.df_du = hbm_derivatives('nl','u',States,hbm,problem);
                
                if isfield(problem,'jacobU')
                    Ju = feval(problem.jacobU,States.df_du,hbm,problem);
                else
                    Ju = sum(hbm.nonlin.hbm.Ju.*States.df_du(ijacobx,ijacobu,:),3);
                end
            end
            varargout{o} = Ju;
        case 'jacobXdot' %df_dX = w*Jxdot*X
            if ~hbm.dependence.xdot
                Jxdot = repmat({zeros(NDofTot,NDofTot)},1,2);
            else
                States.df_dxdot = hbm_derivatives('nl','xdot',States,hbm,problem);

                if isfield(problem,'jacobXdot')
                    Jxdot = feval(problem.jacobXdot,States,hbm,problem);
                else
                    for n = 1:2
                        Jxdot{n} = sum(hbm.nonlin.hbm.Jxdot{n}.*States.df_dxdot(ijacobx,ijacobx,:),3);
                    end
                end
            end
            varargout{o} = Jxdot;
        case 'jacobUdot' %df_dU = w*Judot*U
            if ~hbm.dependence.udot
                Judot = repmat({zeros(NDofTot,NInputTot)},1,2);
            else
                States.df_dudot = hbm_derivatives('nl','udot',States,hbm,problem);
                
                if isfield(problem,'jacobUdot')
                    Judot = feval(problem.jacobUdot,States,hbm,problem);
                else
                    for n = 1:2
                        Judot{n} = sum(hbm.nonlin.hbm.Judot{n}.*States.df_dudot(ijacobx,ijacobu,:),3);
                    end
                end
            end
            varargout{o} = Judot;
        case 'jacobXddot' %df_dX = w^2*Jxddot*X
            if ~hbm.dependence.xddot
                Jxddot = repmat({zeros(NDofTot,NDofTot)},1,3);
            else
                States.df_dxddot = hbm_derivatives('nl','xddot',States,hbm,problem);
                
                if isfield(problem,'jacobXdot')
                    Jxddot = feval(problem.jacobXdot,States,hbm,problem);
                else
                    for n = 1:3
                        Jxddot{n} = sum(hbm.nonlin.hbm.Jxddot{n}.*States.df_dxddot(ijacobx,ijacobx,:),3);
                    end
                end
            end
            varargout{o} = Jxddot;
        case 'jacobUddot' %df_dU = w^2*Juddot*U
            if ~hbm.dependence.uddot
                Juddot = repmat({zeros(NDofTot,NInputTot)},1,3);
            else
                States.df_duddot = hbm_derivatives('nl','uddot',States,hbm,problem);
                
                if isfield(problem,'jacobUddot')
                    Juddot = feval(problem.jacobUddot,States,hbm,problem);
                else
                    for n = 1:3
                        Juddot{n} = sum(hbm.nonlin.hbm.Juddot{n}.*States.df_duddot(ijacobx,ijacobu,:),3);
                    end
                end
            end
            varargout{o} = Juddot;
        case 'floquet1xdot'
           if ~hbm.dependence.xdot
                D1 = zeros(NDofTot,NDofTot);
            else
                States.df_dxdot = hbm_derivatives('nl','xdot',States,hbm,problem);
                
                if isfield(problem,'floquet1xdot')
                    D1 = feval(problem.floquet1xdot,States,hbm,problem);
                else
                    D1 = sum(hbm.nonlin.hbm.Jx.*States.df_dxdot(ijacobx,ijacobx,:),3);
                end
            end
            varargout{o} = D1;
        case 'floquet1xddot'
           if ~hbm.dependence.xddot
                D1dd = repmat({zeros(NDofTot,NDofTot)},1,2);
            else
                States.df_dxddot = hbm_derivatives('nl','xddot',States,hbm,problem);
                
                if isfield(problem,'floquet1xddot')
                    D1dd = feval(problem.floquet1xddot,States,hbm,problem);
                else
                    for n = 1:2
                        D1dd{n} = sum(hbm.nonlin.hbm.Jxdot{n}.*States.df_dxddot(ijacobx,ijacobx,:),3);
                    end
                end
            end
            varargout{o} = D1dd;
       case 'floquet2'
           if ~hbm.dependence.xddot
               D2 = zeros(NDofTot);
           else
               States.df_dxddot = hbm_derivatives('nl','xddot',States,hbm,problem);
               
               if isfield(problem,'floquet2')
                   D2 = feval(problem.floquet2,States,hbm,problem);
               else
                   D2 = sum(hbm.nonlin.hbm.Jx.*States.df_dxddot(ijacobx,ijacobx,:),3);
               end
            end
            varargout{o} = D2;
    end
end