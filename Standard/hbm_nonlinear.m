function varargout = hbm_nonlinear(command,hbm,problem,w0,Xp,Up)
NInput = problem.NInput;
NDof = problem.NDof;

NFreq = hbm.harm.NFreq;
NComp = hbm.harm.NComp;
Nfft  = hbm.harm.Nfft(1);
kHarm = hbm.harm.kHarm(:,1);
wBase = hbm.harm.rFreqBase(1)*w0;

iRetain = hbm.harm.iRetain;
NRetain = hbm.harm.NRetain;

%unpack the inputs
w = kHarm*wBase;
States.t = (0:Nfft-1)/Nfft*2*pi/wBase;

%work out the time domain
X = unpackdof(Xp,NFreq-1,NDof,iRetain);
U = unpackdof(Up,NFreq-1,NInput);

%compute the fourier coefficients of the derivatives
Wx = repmat(1i*w,1,size(X,2));
Xdot  = X.*Wx;
Xddot = Xdot.*Wx;

%precompute the external inputs
Wu = repmat(1i*w,1,size(U,2));
Udot  = U.*Wu;
Uddot = Udot.*Wu;

States.w0 = w0;
States.wBase = wBase;
switch hbm.options.aft_method
    case 'fft'
        %create the time series from the fourier series
        States.x     = freq2time(X    ,NFreq-1,Nfft).';
        States.xdot  = freq2time(Xdot ,NFreq-1,Nfft).';
        States.xddot = freq2time(Xddot,NFreq-1,Nfft).';
        
        %create the vector of inputs
        States.u     = freq2time(U    ,NFreq-1,Nfft).';
        States.udot  = freq2time(Udot ,NFreq-1,Nfft).';
        States.uddot = freq2time(Uddot,NFreq-1,Nfft).';
        
    case 'mat'
        %create the time series from the fourier series
        States.x     = real(hbm.nonlin.IFFT*X).';
        States.xdot  = real(hbm.nonlin.IFFT*Xdot).';
        States.xddot = real(hbm.nonlin.IFFT*Xddot).';
        
        %create the vector of inputs
        States.u     = real(hbm.nonlin.IFFT*U).';
        States.udot  = real(hbm.nonlin.IFFT*Udot).';
        States.uddot = real(hbm.nonlin.IFFT*Uddot).';
end

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

Fp = packdof(F,iRetain);
  
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
                Dw = zeros(NRetain,1);
            else
                dfhbm_dw = hbm_derivatives('nl' ,'w',States,hbm,problem);
                Dw = packdof(hbm.nonlin.FFT*dfhbm_dw{1},iRetain); 
            end
            varargout{o} = Dw;
        case 'jacobX'  %df_dX = Jx*X
            if ~hbm.dependence.x
                Jx = zeros(NRetain);
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
                            
                            Jx = zeros(NComp*NDof);
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
                            Jx = Jx(iRetain,iRetain);
                    end
                end
            end
            varargout{o} = Jx;
        case 'jacobU' %df_dU = Ju*U
            if ~hbm.dependence.u
                Ju = zeros(NRetain,NComp*NInput);
            else
                States.df_du = hbm_derivatives('nl' ,'u',States,hbm,problem);
                
                if isfield(problem,'jacobU')
                    Ju = feval(problem.jacobU,States,problem,w0);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Ju = sum(hbm.nonlin.hbm.Ju.*States.df_du(ijacobx,ijacobu,:),3);
                        case 'sum'
                            theta = permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]);
                            
                            Ju = zeros(NComp*NInput);
                            Ju(1:NInput,1:NInput) = mean(States.df_du,3);
                            for l = 1:(NFreq-1)
                                cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                Ju(1:NInput,(2*l-1)*NInput+(1:NInput)) = sum( (States.df_du.*cl)/Nfft,3);
                                Ju(1:NInput,(2*l-0)*NInput+(1:NInput)) = sum(-(States.df_du.*sl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                ck = cos(kHarm(k+1,1)*theta); sk = sin(kHarm(k+1,1)*theta);
                                Ju((2*k-1)*NInput+(1:NInput),1:NInput) =  2*sum(States.df_du.*ck/Nfft,3);
                                Ju((2*k-0)*NInput+(1:NInput),1:NInput) = -2*sum(States.df_du.*sk/Nfft,3);
                                for l = 1:(NFreq-1)
                                    cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                    Ju((2*k-1)*NInput+(1:NInput),(2*l-1)*NInput+(1:NInput)) =  2*sum( (States.df_du.*cl).*ck/Nfft,3);
                                    Ju((2*k-1)*NInput+(1:NInput),(2*l-0)*NInput+(1:NInput)) =  2*sum(-(States.df_du.*sl).*ck/Nfft,3);
                                    Ju((2*k-0)*NInput+(1:NInput),(2*l-1)*NInput+(1:NInput)) = -2*sum( (States.df_du.*cl).*sk/Nfft,3);
                                    Ju((2*k-0)*NInput+(1:NInput),(2*l-0)*NInput+(1:NInput)) = -2*sum(-(States.df_du.*sl).*sk/Nfft,3);
                                end
                            end
                            Ju = Ju(iRetain,:);
                    end
                end
            end
            varargout{o} = Ju;
        case 'jacobXdot' %df_dX = w*Jxdot*X
             if ~hbm.dependence.xdot
                 Jxdot = zeros(NRetain);
             else
                 States.df_dxdot = hbm_derivatives('nl' ,'xdot',States,hbm,problem);
                 
                 if isfield(problem,'jacobXdot')
                     Jxdot = feval(problem.jacobXdot,States,hbm,problem);
                 else
                     switch hbm.options.jacob_method
                         case 'mat'
                             Jxdot = sum(hbm.nonlin.hbm.Jxdot{1}.*States.df_dxdot(ijacobx,ijacobx,:),3);
                         case 'sum'
                             theta = permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]);
                             
                             Jxdot = zeros(NComp*NDof);
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
                     Jxdot = Jxdot(iRetain,iRetain);
                 end
             end
            varargout{o} = Jxdot;
        case 'jacobUdot' %df_dU = w*Judot*U
            if ~hbm.dependence.udot
                Judot = zeros(NRetain,NComp*NInput);
            else
                States.df_dudot = hbm_derivatives('nl' ,'udot',States,hbm,problem);
                
                if isfield(problem,'jacobUdot')
                    Judot = feval(problem.jacobUdot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Judot = sum(hbm.nonlin.hbm.Judot{1}.*States.df_dudot(ijacobx,ijacobu,:),3);
                        case 'sum'
                            theta = permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]);
                            
                            Judot = zeros(NComp*NDof,NComp*NInput);
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
                    Judot = Judot(iRetain,:);
                end
            end
            varargout{o} = Judot;
         case 'jacobXddot' %df_dX = w^2*Jxddot*X
            if ~hbm.dependence.xddot
                Jxddot = zeros(NRetain);
            else
                States.df_dxddot = hbm_derivatives('nl' ,'xddot',States,hbm,problem);
                
                if isfield(problem,'jacobXddot')
                    Jxddot = feval(problem.jacobXddot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Jxddot = sum(hbm.nonlin.hbm.Jxddot{1}.*States.df_dxddot(ijacobx,ijacobx,:),3);
                        case 'sum'
                            theta = permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]);
                            
                            Jxddot = zeros(NDof*NComp);
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
                    Jxddot = Jxddot(iRetain,iRetain);
                end
            end
            varargout{o} = Jxddot;
        case 'jacobUddot' %df_dU = w^2*Juddot*U
            if ~hbm.dependence.uddot
                Juddot = zeros(NRetain,NComp*NInput);
            else
                States.df_duddot = hbm_derivatives('nl' ,'uddot',States,hbm,problem);
                
                if isfield(problem,'jacobUddot')
                    Juddot = feval(problem.jacobUddot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Juddot = sum(hbm.nonlin.hbm.Juddot{1}.*States.df_duddot(ijacobx,ijacobu,:),3);
                        case 'sum'
                            theta = permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]);
                            
                            Juddot = zeros(NDof*NComp);
                            for l = 1:(NFreq-1)
                                [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                Juddot(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum((-kHarm(l+1,1)^2*States.df_duddot.*cl)/Nfft,3);
                                Juddot(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum(( kHarm(l+1,1)^2*States.df_duddot.*sl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                [ck, sk] = cosAndSin(kHarm(k+1,1)*theta);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                    Juddot((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum((- kHarm(l+1,1)^2*States.df_duddot.*cl).*ck/Nfft,3);
                                    Juddot((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((  kHarm(l+1,1)^2*States.df_duddot.*sl).*ck/Nfft,3);
                                    Juddot((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum((- kHarm(l+1,1)^2*States.df_duddot.*cl).*sk/Nfft,3);
                                    Juddot((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((  kHarm(l+1,1)^2*States.df_duddot.*sl).*sk/Nfft,3);
                                end
                            end
                    end
                    Juddot = Juddot(iRetain,:);
                end
            end
            varargout{o} = Juddot;
            
       case 'floquet1xdot'
           if ~hbm.dependence.xdot
                D1 = zeros(NRetain);
            else
                States.df_dxdot = hbm_derivatives('nl' ,'xdot',States,hbm,problem);
                
                if isfield(problem,'floquet1xdot')
                    D1 = feval(problem.floquet1xdot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            D1 = sum(hbm.nonlin.hbm.Jx.*States.df_dxdot(ijacobx,ijacobx,:),3);
                        case 'sum'
                            theta = permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]);
                            
                            D1 = zeros(NDof*NComp);
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
                    D1 = D1(iRetain,iRetain);
                end
            end
            varargout{o} = D1;
        case 'floquet1xddot'
           if ~hbm.dependence.xddot
                D1 = zeros(NRetain);
            else
                States.df_dxddot = hbm_derivatives('nl' ,'xddot',States,hbm,problem);
                
                if isfield(problem,'floquet1xddot')
                    D1 = feval(problem.floquet1xddot,States,hbm,problem);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            D1 = sum(hbm.nonlin.hbm.Jxdot{1}.*States.df_dxddot(ijacobx,ijacobx,:),3);
                        case 'sum'
                            theta = repmat(permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]),NDof,NDof,1);
                            
                            D1 = zeros(NDof*NComp);
                            for l = 1:(NFreq-1)
                                [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                D1(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum((-kHarm(l+1,1)*States.df_dxddot.*sl)/Nfft,3);
                                D1(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum((-kHarm(l+1,1)*States.df_dxddot.*cl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                [ck, sk] = cosAndSin(kHarm(k+1,1)*theta);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                    D1((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum((-kHarm(l+1,1)*States.df_dxddot.*sl).*ck/Nfft,3);
                                    D1((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((-kHarm(l+1,1)*States.df_dxddot.*cl).*ck/Nfft,3);
                                    D1((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum((-kHarm(l+1,1)*States.df_dxddot.*sl).*sk/Nfft,3);
                                    D1((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((-kHarm(l+1,1)*States.df_dxddot.*cl).*sk/Nfft,3);
                                end
                            end
                    end
                    D1 = D1(iRetain,iRetain);
                end
            end
            varargout{o} = D1;
            
       case 'floquet2'
           if ~hbm.dependence.xddot
                D2 = zeros(NRetain);
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
                            
                            D2 = zeros(NDof*NComp);
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
                    D2 = D2(iRetain,iRetain);
                end
            end
            varargout{o} = D2;
           
    end
end