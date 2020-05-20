function varargout = hbm_nonlinear(command,hbm,problem,w0,Xp,Up)
NInput = problem.NInput;
NDof = problem.NDof;

ii = find(hbm.harm.NHarm ~= 0);

NFreq = hbm.harm.NFreq;
NComp = hbm.harm.NComp;
NDofTot = NComp*NDof;
NInputTot = NComp*NInput;

%work out the time domain
X = unpackdof(Xp,NFreq-1,NDof);
U = unpackdof(Up,NFreq-1,NInput);

%get the time series
States = hbm_states(w0,X,U,hbm);

%push through the nl system
States.f = feval(problem.model,'nl' ,States,hbm,problem);

%finally convert into a fourier series
switch hbm.options.aft_method
    case 'fft'
        F = time2freq(States.f.',NFreq-1,Nfft);
    case 'mat'
        %finally convert into a fourier series
        F = hbm.nonlin.FFT*States.f.';
end

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
                Dw = zeros(NDofTot,1);
            else
                dfhbm_dw = hbm_derivatives('nl' ,'w',States,hbm,problem);
                if isfield(problem,'derivW')
                    Dw = feval(problem.derivW,dfhbm_dw,hbm,problem);
                else
                    Dw = packdof(hbm.nonlin.FFT*dfhbm_dw{ii});
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
                    switch hbm.options.jacob_method
                        case 'mat'
                            Jx = sum(hbm.nonlin.hbm.Jx.*States.df_dx(ijacobx,ijacobx,:),3);
                        case 'sum'
                            theta = permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]);
                            
                            Jx = zeros(NDofTot,NDofTot);
                            Jx(1:NDof,1:NDof) = mean(States.df_dx,3);
                            for l = 1:(NFreq-1)
                                cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                Jx(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum( (States.df_dx.*cl)/Nfft,3);
                                Jx(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum(-(States.df_dx.*sl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                ck = cos(kHarm(k+1,1)*theta); sk = sin(kHarm(k+1,1)*theta);
                                Jx((2*k-1)*NDof+(1:NDof),1:NDof) =  2*sum(States.df_dx.*ck/Nfft,3);
                                Jx((2*k-0)*NDof+(1:NDof),1:NDof) = -2*sum(States.df_dx.*sk/Nfft,3);
                                for l = 1:(NFreq-1)
                                    cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                    Jx((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum( (States.df_dx.*cl).*ck/Nfft,3);
                                    Jx((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum(-(States.df_dx.*sl).*ck/Nfft,3);
                                    Jx((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum( (States.df_dx.*cl).*sk/Nfft,3);
                                    Jx((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum(-(States.df_dx.*sl).*sk/Nfft,3);
                                end
                            end
                    end
                end
            end
            varargout{o} = Jx;
        case 'jacobU' %df_dU = Ju*U
            if ~hbm.dependence.u
                Ju = zeros(NDofTot,NInputTot);
            else
                States.df_du = hbm_derivatives('nl' ,'u',States,hbm,problem);
                
                if isfield(problem,'jacobU')
                    Ju = feval(problem.jacobU,States,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Ju = sum(hbm.nonlin.hbm.Ju.*States.df_du(ijacobx,ijacobu,:),3);
                        case 'sum'
                            theta = permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]);
                            
                            Ju = zeros(NDofTot,NInputTot);
                            Ju(1:NDof,1:NInput) = mean(States.df_du,3);
                            for l = 1:(NFreq-1)
                                cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                Ju(1:NDof,(2*l-1)*NInput+(1:NInput)) = sum( (States.df_du.*cl)/Nfft,3);
                                Ju(1:NDof,(2*l-0)*NInput+(1:NInput)) = sum(-(States.df_du.*sl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                ck = cos(kHarm(k+1,1)*theta); sk = sin(kHarm(k+1,1)*theta);
                                Ju((2*k-1)*NDof+(1:NDof),1:NInput) =  2*sum(States.df_du.*ck/Nfft,3);
                                Ju((2*k-0)*NDof+(1:NDof),1:NInput) = -2*sum(States.df_du.*sk/Nfft,3);
                                for l = 1:(NFreq-1)
                                    cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                    Ju((2*k-1)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) =  2*sum( (States.df_du.*cl).*ck/Nfft,3);
                                    Ju((2*k-1)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) =  2*sum(-(States.df_du.*sl).*ck/Nfft,3);
                                    Ju((2*k-0)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) = -2*sum( (States.df_du.*cl).*sk/Nfft,3);
                                    Ju((2*k-0)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) = -2*sum(-(States.df_du.*sl).*sk/Nfft,3);
                                end
                            end
                    end
                end
            end
            varargout{o} = Ju;
        case 'jacobXdot' %df_dX = w*Jxdot*X
             if ~hbm.dependence.xdot
                 Jxdot = zeros(DofTot,NDofTot);
             else
                 States.df_dxdot = hbm_derivatives('nl' ,'xdot',States,hbm,problem);
                 
                 if isfield(problem,'jacobXdot')
                     Jxdot = feval(problem.jacobXdot,States,hbm,problem);
                 else
                     switch hbm.options.jacob_method
                         case 'mat'
                             Jxdot = sum(hbm.nonlin.hbm.Jxdot{ii}.*States.df_dxdot(ijacobx,ijacobx,:),3);
                         case 'sum'
                             theta = permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]);
                             
                             Jxdot = zeros(NDofTot,NDofTot);
                             for l = 1:(NFreq-1)
                                 cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                 Jxdot(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum((-kHarm(l+1,1)*States.df_dxdot.*sl)/Nfft,3);
                                 Jxdot(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum((-kHarm(l+1,1)*States.df_dxdot.*cl)/Nfft,3);
                             end
                             for k = 1:(NFreq-1)
                                 ck = cos(kHarm(k+1,1)*theta); sk = sin(kHarm(k+1,1)*theta);
                                 for l = 1:(NFreq-1)
                                     cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                     Jxdot((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum((- kHarm(l+1,1)*States.df_dxdot.*sl).*ck/Nfft,3);
                                     Jxdot((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((- kHarm(l+1,1)*States.df_dxdot.*cl).*ck/Nfft,3);
                                     Jxdot((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum((- kHarm(l+1,1)*States.df_dxdot.*sl).*sk/Nfft,3);
                                     Jxdot((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((- kHarm(l+1,1)*States.df_dxdot.*cl).*sk/Nfft,3);
                                 end
                             end
                     end
                 end
             end
            varargout{o} = Jxdot;
        case 'jacobUdot' %df_dU = w*Judot*U
            if ~hbm.dependence.udot
                Judot = zeros(NDofTot,NInputTot);
            else
                States.df_dudot = hbm_derivatives('nl' ,'udot',States,hbm,problem);
                
                if isfield(problem,'jacobUdot')
                    Judot = feval(problem.jacobUdot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Judot = sum(hbm.nonlin.hbm.Judot{ii}.*States.df_dudot(ijacobx,ijacobu,:),3);
                        case 'sum'
                            theta = permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]);
                            
                            Judot = zeros(NDofTot,NInputTot);
                            for l = 1:(NFreq-1)
                                cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                Judot(1:NDof,(2*l-1)*NInput+(1:NInput)) =  sum((-kHarm(l+1,1)*States.df_dudot.*sl)/Nfft,3);
                                Judot(1:NDof,(2*l-0)*NInput+(1:NInput)) =  sum((-kHarm(l+1,1)*States.df_dudot.*cl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                ck = cos(kHarm(k+1,1)*theta); sk = sin(kHarm(k+1,1)*theta);
                                for l = 1:(NFreq-1)
                                    cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                    Judot((2*k-1)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) =  2*sum((- kHarm(l+1,1)*States.df_dudot.*sl).*ck/Nfft,3);
                                    Judot((2*k-1)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) =  2*sum((- kHarm(l+1,1)*States.df_dudot.*cl).*ck/Nfft,3);
                                    Judot((2*k-0)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) = -2*sum((- kHarm(l+1,1)*States.df_dudot.*sl).*sk/Nfft,3);
                                    Judot((2*k-0)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) = -2*sum((- kHarm(l+1,1)*States.df_dudot.*cl).*sk/Nfft,3);
                                end
                            end
                    end
                end
            end
            varargout{o} = Judot;
         case 'jacobXddot' %df_dX = w^2*Jxddot*X
            if ~hbm.dependence.xddot
                Jxddot = zeros(NDofTot,NDofTot);
            else
                States.df_dxddot = hbm_derivatives('nl' ,'xddot',States,hbm,problem);
                
                if isfield(problem,'jacobXddot')
                    Jxddot = feval(problem.jacobXddot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Jxddot = sum(hbm.nonlin.hbm.Jxddot{ii}.*States.df_dxddot(ijacobx,ijacobx,:),3);
                        case 'sum'
                            theta = permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]);
                            
                            Jxddot = zeros(NDofTot,NDofTot);
                            for l = 1:(NFreq-1)
                                [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                Jxddot(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum((-kHarm(l+1,1)^2*States.df_dxddot.*cl)/Nfft,3);
                                Jxddot(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum(( kHarm(l+1,1)^2*States.df_dxddot.*sl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                [ck, sk] = cosAndSin(kHarm(k+1,1)*theta);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                    Jxddot((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum((- kHarm(l+1,1)^2*States.df_dxddot.*cl).*ck/Nfft,3);
                                    Jxddot((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((  kHarm(l+1,1)^2*States.df_dxddot.*sl).*ck/Nfft,3);
                                    Jxddot((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum((- kHarm(l+1,1)^2*States.df_dxddot.*cl).*sk/Nfft,3);
                                    Jxddot((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((  kHarm(l+1,1)^2*States.df_dxddot.*sl).*sk/Nfft,3);
                                end
                            end
                    end
                end
            end
            varargout{o} = Jxddot;
        case 'jacobUddot' %df_dU = w^2*Juddot*U
            if ~hbm.dependence.uddot
                Juddot = zeros(NDofTot,NInputTot);
            else
                States.df_duddot = hbm_derivatives('nl' ,'uddot',States,hbm,problem);
                
                if isfield(problem,'jacobUddot')
                    Juddot = feval(problem.jacobUddot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Juddot = sum(hbm.nonlin.hbm.Juddot{ii}.*States.df_duddot(ijacobx,ijacobu,:),3);
                        case 'sum'
                            theta = permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]);
                            
                            Juddot = zeros(NDofTot,NInputTot);
                            for l = 1:(NFreq-1)
                                [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                Juddot(1:NDof,(2*l-1)*NInput+(1:NInput)) =  sum((-kHarm(l+1,1)^2*States.df_duddot.*cl)/Nfft,3);
                                Juddot(1:NDof,(2*l-0)*NInput+(1:NInput)) =  sum(( kHarm(l+1,1)^2*States.df_duddot.*sl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                [ck, sk] = cosAndSin(kHarm(k+1,1)*theta);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                    Juddot((2*k-1)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) =  2*sum((- kHarm(l+1,1)^2*States.df_duddot.*cl).*ck/Nfft,3);
                                    Juddot((2*k-1)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) =  2*sum((  kHarm(l+1,1)^2*States.df_duddot.*sl).*ck/Nfft,3);
                                    Juddot((2*k-0)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) = -2*sum((- kHarm(l+1,1)^2*States.df_duddot.*cl).*sk/Nfft,3);
                                    Juddot((2*k-0)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) = -2*sum((  kHarm(l+1,1)^2*States.df_duddot.*sl).*sk/Nfft,3);
                                end
                            end
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
                    switch hbm.options.jacob_method
                        case 'mat'
                            D1 = sum(hbm.nonlin.hbm.Jx.*States.df_dxdot(ijacobx,ijacobx,:),3);
                        case 'sum'
                            theta = permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]);
                            
                            D1 = zeros(NDofTot,NDofTot);
                            D1(1:NDof,1:NDof) = mean(States.df_dxdot,3);
                            for l = 1:(NFreq-1)
                                [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                D1(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum(( States.df_dxdot.*cl)/Nfft,3);
                                D1(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum((-States.df_dxdot.*sl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                [ck, sk] = cosAndSin(kHarm(k+1,1)*theta);
                                D1((2*k-1)*NDof+(1:NDof),1:NDof) =  2*sum(States.df_dxdot.*ck/Nfft,3);
                                D1((2*k-0)*NDof+(1:NDof),1:NDof) = -2*sum(States.df_dxdot.*sk/Nfft,3);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                    D1((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum(( States.df_dxdot.*cl).*ck/Nfft,3);
                                    D1((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((-States.df_dxdot.*sl).*ck/Nfft,3);
                                    D1((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum(( States.df_dxdot.*cl).*sk/Nfft,3);
                                    D1((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((-States.df_dxdot.*sl).*sk/Nfft,3);
                                end
                            end
                    end
                end
            end
            varargout{o} = D1;
        case 'floquet1xddot'
           if ~hbm.dependence.xddot
                D1dd = zeros(NDofTot,NDofTot);
            else
                States.df_dxddot = hbm_derivatives('nl' ,'xddot',States,hbm,problem);
                
                if isfield(problem,'floquet1xddot')
                    D1dd = feval(problem.floquet1xddot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            D1dd = sum(hbm.nonlin.hbm.Jxdot{ii}.*States.df_dxddot(ijacobx,ijacobx,:),3);
                        case 'sum'
                            theta = repmat(permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]),NDof,NDof,1);
                            
                            D1dd = zeros(NDofTot,NDofTot);
                            for l = 1:(NFreq-1)
                                [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                D1dd(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum((-kHarm(l+1,1)*States.df_dxddot.*sl)/Nfft,3);
                                D1dd(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum((-kHarm(l+1,1)*States.df_dxddot.*cl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                [ck, sk] = cosAndSin(kHarm(k+1,1)*theta);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                    D1dd((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum((-kHarm(l+1,1)*States.df_dxddot.*sl).*ck/Nfft,3);
                                    D1dd((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((-kHarm(l+1,1)*States.df_dxddot.*cl).*ck/Nfft,3);
                                    D1dd((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum((-kHarm(l+1,1)*States.df_dxddot.*sl).*sk/Nfft,3);
                                    D1dd((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((-kHarm(l+1,1)*States.df_dxddot.*cl).*sk/Nfft,3);
                                end
                            end
                    end
                end
            end
            varargout{o} = D1dd; 
       case 'floquet2'
           if ~hbm.dependence.xddot
                D2 = zeros(NDofTot,NDofTot);
            else
                States.df_dxddot = hbm_derivatives('nl' ,'xddot',States,hbm,problem);
                
                if isfield(problem,'floquet2')
                    D2 = feval(problem.floquet2,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            D2 = sum(hbm.nonlin.hbm.Jx.*States.df_dxddot(ijacobx,ijacobx,:),3);
                        case 'sum'
                            theta = permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]);
                            
                            D2 = zeros(NDofTot,NDofTot);
                            D2(1:NDof,1:NDof) = mean(States.df_dxddot,3);
                            for l = 1:(NFreq-1)
                                [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                D2(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum(( States.df_dxddot.*cl)/Nfft,3);
                                D2(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum((-States.df_dxddot.*sl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                [ck, sk] = cosAndSin(kHarm(k+1,1)*theta);
                                D2((2*k-1)*NDof+(1:NDof),1:NDof) =  2*sum(States.df_dxddot.*ck/Nfft,3);
                                D2((2*k-0)*NDof+(1:NDof),1:NDof) = -2*sum(States.df_dxddot.*sk/Nfft,3);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                    D2((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum(( States.df_dxddot.*cl).*ck/Nfft,3);
                                    D2((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((-States.df_dxddot.*sl).*ck/Nfft,3);
                                    D2((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum(( States.df_dxddot.*cl).*sk/Nfft,3);
                                    D2((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((-States.df_dxddot.*sl).*sk/Nfft,3);
                                end
                            end
                    end
                end
            end
            varargout{o} = D2;
    end
end