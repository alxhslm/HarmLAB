function varargout = hbm_nonlinear(command,hbm,problem,w0,Xp,Up)
NInput = problem.NInput;
NDof = problem.NDof;
NAlg   = problem.NAlg;

NFreq = hbm.harm.NFreq;
NComp = hbm.harm.NComp;
Nfft  = hbm.harm.Nfft(1);
kHarm = hbm.harm.kHarm(:,1);
rBase = hbm.harm.rFreqBase(1);

iRetain = hbm.harm.iRetain;
NRetain = hbm.harm.NRetain;

%unpack the inputs
w = kHarm*rBase*w0;
t = (0:Nfft-1)'/Nfft*2*pi/(rBase*w0);

%work out the time domain
X = unpackdof(Xp,NFreq-1,NDof,iRetain);
U = unpackdof(Up,NFreq-1,NInput);
xalg = reshape(Xp(NRetain+1:end),NAlg,prod(Nfft))';

%compute the fourier coefficients of the derivatives
Wx = repmat(1i*w,1,size(X,2));
Xdot  = X.*Wx;
Xddot = Xdot.*Wx;

%precompute the external inputs
Wu = repmat(1i*w,1,size(U,2));
Udot  = U.*Wu;
Uddot = Udot.*Wu;

switch hbm.options.aft_method
    case 'fft'
        %create the time series from the fourier series
        x     = freq2time(X    ,NFreq-1,Nfft);
        xdot  = freq2time(Xdot ,NFreq-1,Nfft);
        xddot = freq2time(Xddot,NFreq-1,Nfft);
        
        %create the vector of inputs
        u     = freq2time(U    ,NFreq-1,Nfft);
        udot  = freq2time(Udot ,NFreq-1,Nfft);
        uddot = freq2time(Uddot,NFreq-1,Nfft);
        
    case 'mat'
        %create the time series from the fourier series
        x     = real(hbm.nonlin.IFFT*X);
        xdot  = real(hbm.nonlin.IFFT*Xdot);
        xddot = real(hbm.nonlin.IFFT*Xddot);
        
        %create the vector of inputs
        u     = real(hbm.nonlin.IFFT*U);
        udot  = real(hbm.nonlin.IFFT*Udot);
        uddot = real(hbm.nonlin.IFFT*Uddot);
end

%push through the nl system
f_hbm = feval(problem.model,'nl' ,t',x.',xdot.',xddot.',u.',udot.',uddot.',xalg.',hbm,problem,w0).';
f_alg = feval(problem.model,'alg',t',x.',xdot.',xddot.',u.',udot.',uddot.',xalg.',hbm,problem,w0).';

%finally convert into a fourier series
switch hbm.options.aft_method
    case 'fft'
        F = time2freq(f_hbm,NFreq-1,Nfft);
    case 'mat'
        %finally convert into a fourier series
        F = hbm.nonlin.FFT*f_hbm;
end

Fp = [packdof(F,iRetain);
      reshape(f_alg',[],1)];
  
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
                dfhbm_dw = hbm_derivatives('nl' ,'w',t,x,xdot,xddot,u,udot,uddot,xalg,f_hbm,hbm,problem,w0);
                dfalg_dw = hbm_derivatives('alg','w',t,x,xdot,xddot,u,udot,uddot,xalg,f_alg,hbm,problem,w0);
                Dw = packdof(hbm.nonlin.FFT*dfhbm_dw{1},iRetain); 
                if problem.NAlg > 0
                    Dw = [Dw; catmat(dfalg_dw{1},1)];
                end
            end
            varargout{o} = Dw;
        case 'jacobX'  %df_dX = Jx*X
            if ~hbm.dependence.x
                Jx = zeros(NRetain);
            else
                dfhbm_dxhbm = hbm_derivatives('nl','x'  ,t,x,xdot,xddot,u,udot,uddot,xalg,f_hbm,hbm,problem,w0);
                dfhbm_dxalg = hbm_derivatives('nl','alg',t,x,xdot,xddot,u,udot,uddot,xalg,f_hbm,hbm,problem,w0);
                
                dfalg_dxhbm = hbm_derivatives('alg','x'  ,t,x,xdot,xddot,u,udot,uddot,xalg,f_alg,hbm,problem,w0);
                dfalg_dxalg = hbm_derivatives('alg','alg',t,x,xdot,xddot,u,udot,uddot,xalg,f_alg,hbm,problem,w0);
                
                if isfield(problem,'jacobX')
                    Jx = feval(problem.jacobX,dfhbm_dxhbm,dfhbm_dxalg,dfalg_dxhbm,dfalg_dxalg,hbm,problem,w0);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Jxx = sum(hbm.nonlin.hbm.Jx.*dfhbm_dxhbm(ijacobx,ijacobx,:),3);
                            Jxa = catmat(hbm.nonlin.alg.Jf.*dfhbm_dxalg(ijacobx,:,:),2);
                            Jax = catmat(hbm.nonlin.alg.Jx.*dfalg_dxhbm(:,ijacobx,:),1);
                            Jaa = blkmat(dfalg_dxalg);
                            Jx = [ Jxx   Jxa;
                                Jax   Jaa];
                        case 'sum'
                            theta = repmat(permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]),NDof,NDof,1);
                            Jx = zeros(NComp*NDof);
                            Jx(1:NDof,1:NDof) = mean(dfhbm_dxhbm,3);
                            for l = 1:(NFreq-1)
                                cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                Jx(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum( (dfhbm_dxhbm.*cl)/Nfft,3);
                                Jx(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum(-(dfhbm_dxhbm.*sl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                ck = cos(kHarm(k+1,1)*theta); sk = sin(kHarm(k+1,1)*theta);
                                Jx((2*k-1)*NDof+(1:NDof),1:NDof) =  2*sum(dfhbm_dxhbm.*ck/Nfft,3);
                                Jx((2*k-0)*NDof+(1:NDof),1:NDof) = -2*sum(dfhbm_dxhbm.*sk/Nfft,3);
                                for l = 1:(NFreq-1)
                                    cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                    Jx((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum( (dfhbm_dxhbm.*cl).*ck/Nfft,3);
                                    Jx((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum(-(dfhbm_dxhbm.*sl).*ck/Nfft,3);
                                    Jx((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum( (dfhbm_dxhbm.*cl).*sk/Nfft,3);
                                    Jx((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum(-(dfhbm_dxhbm.*sl).*sk/Nfft,3);
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
                dfhbm_du = hbm_derivatives('nl' ,'u',t,x,xdot,xddot,u,udot,uddot,xalg,f_hbm,hbm,problem,w0);
                dfalg_du = hbm_derivatives('alg','u',t,x,xdot,xddot,u,udot,uddot,xalg,f_alg,hbm,problem,w0);
                
                if isfield(problem,'jacobU')
                    Ju = feval(problem.jacobU,dfhbm_du,dfalg_du,hbm,problem,w0);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Jxu = sum(hbm.nonlin.hbm.Ju.*dfhbm_du(ijacobx,ijacobu,:),3);
                            Jau = catmat(hbm.nonlin.alg.Ju.*dfalg_du(:,ijacobu,:),1);
                            Ju = [Jxu; Jau];
                        case 'sum'
                            theta = repmat(permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]),NDof,NDof,1);
                            Ju = zeros(NComp*NInput);
                            Ju(1:NInput,1:NInput) = mean(dfhbm_du,3);
                            for l = 1:(NFreq-1)
                                cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                Ju(1:NInput,(2*l-1)*NInput+(1:NInput)) = sum( (dfhbm_du.*cl)/Nfft,3);
                                Ju(1:NInput,(2*l-0)*NInput+(1:NInput)) = sum(-(dfhbm_du.*sl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                ck = cos(kHarm(k+1,1)*theta); sk = sin(kHarm(k+1,1)*theta);
                                Ju((2*k-1)*NInput+(1:NInput),1:NInput) =  2*sum(dfhbm_du.*ck/Nfft,3);
                                Ju((2*k-0)*NInput+(1:NInput),1:NInput) = -2*sum(dfhbm_du.*sk/Nfft,3);
                                for l = 1:(NFreq-1)
                                    cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                    Ju((2*k-1)*NInput+(1:NInput),(2*l-1)*NInput+(1:NInput)) =  2*sum( (dfhbm_du.*cl).*ck/Nfft,3);
                                    Ju((2*k-1)*NInput+(1:NInput),(2*l-0)*NInput+(1:NInput)) =  2*sum(-(dfhbm_du.*sl).*ck/Nfft,3);
                                    Ju((2*k-0)*NInput+(1:NInput),(2*l-1)*NInput+(1:NInput)) = -2*sum( (dfhbm_du.*cl).*sk/Nfft,3);
                                    Ju((2*k-0)*NInput+(1:NInput),(2*l-0)*NInput+(1:NInput)) = -2*sum(-(dfhbm_du.*sl).*sk/Nfft,3);
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
                 dfhbm_dxdot = hbm_derivatives('nl' ,'xdot',t,x,xdot,xddot,u,udot,uddot,xalg,f_hbm,hbm,problem,w0);
                 dfalg_dxdot = hbm_derivatives('alg','xdot',t,x,xdot,xddot,u,udot,uddot,xalg,f_alg,hbm,problem,w0);
                 
                 if isfield(problem,'jacobXdot')
                     Jxdot = feval(problem.jacobXdot,dfhbm_dxdot,dfalg_dxdot,hbm,problem,w0);
                 else
                     switch hbm.options.jacob_method
                         case 'mat'
                             Jxx = sum(hbm.nonlin.hbm.Jxdot{1}.*dfhbm_dxdot(ijacobx,ijacobx,:),3);
                             Jxa = zeros(NComp*NDof,NAlg*prod(Nfft));
                             Jax = catmat(hbm.nonlin.alg.Jxdot{1}.*dfalg_dxdot(:,ijacobx,:),1);
                             Jaa = zeros(NAlg*prod(Nfft));
                             Jxdot = [Jxx Jxa;
                                 Jax Jaa];
                         case 'sum'
                             theta = repmat(permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]),NDof,NDof,1);
                             
                             Jxdot = zeros(NComp*NDof);
                             for l = 1:(NFreq-1)
                                 cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                 Jxdot(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum((-kHarm(l+1,1)*dfhbm_dxdot.*sl)/Nfft,3);
                                 Jxdot(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum((-kHarm(l+1,1)*dfhbm_dxdot.*cl)/Nfft,3);
                             end
                             for k = 1:(NFreq-1)
                                 ck = cos(kHarm(k+1,1)*theta); sk = sin(kHarm(k+1,1)*theta);
                                 for l = 1:(NFreq-1)
                                     cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                     Jxdot((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum((- kHarm(l+1,1)*dfhbm_dxdot.*sl).*ck/Nfft,3);
                                     Jxdot((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((- kHarm(l+1,1)*dfhbm_dxdot.*cl).*ck/Nfft,3);
                                     Jxdot((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum((- kHarm(l+1,1)*dfhbm_dxdot.*sl).*sk/Nfft,3);
                                     Jxdot((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((- kHarm(l+1,1)*dfhbm_dxdot.*cl).*sk/Nfft,3);
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
                dfhbm_dudot = hbm_derivatives('nl' ,'udot',t,x,xdot,xddot,u,udot,uddot,xalg,f_hbm,hbm,problem,w0);
                dfalg_dudot = hbm_derivatives('alg','udot',t,x,xdot,xddot,u,udot,uddot,xalg,f_alg,hbm,problem,w0);
                
                if isfield(problem,'jacobUdot')
                    Judot = feval(problem.jacobUdot,dfhbm_dudot,dfalg_dudot,hbm,problem,w0);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Jxu = sum(hbm.nonlin.hbm.Judot{1}.*dfhbm_dudot(ijacobx,ijacobu,:),3);
                            Jau = catmat(hbm.nonlin.alg.Judot{1}.*dfalg_dudot(:,ijacobu,:),1);
                            Judot = [Jxu; Jau];
                        case 'sum'
                            theta = repmat(permute(2*pi/Nfft*(0:(Nfft-1)),[1 3 2]),NDof,NDof,1);
                            
                            Judot = zeros(NComp*NDof,NComp*NInput);
                            for l = 1:(NFreq-1)
                                cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                Judot(1:NDof,(2*l-1)*NInput+(1:NInput)) =  sum((-kHarm(l+1,1)*dfhbm_dudot.*sl)/Nfft,3);
                                Judot(1:NDof,(2*l-0)*NInput+(1:NInput)) =  sum((-kHarm(l+1,1)*dfhbm_dudot.*cl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                ck = cos(kHarm(k+1,1)*theta); sk = sin(kHarm(k+1,1)*theta);
                                for l = 1:(NFreq-1)
                                    cl = cos(kHarm(l+1,1)*theta); sl = sin(kHarm(l+1,1)*theta);
                                    Judot((2*k-1)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) =  2*sum((- kHarm(l+1,1)*dfhbm_dudot.*sl).*ck/Nfft,3);
                                    Judot((2*k-1)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) =  2*sum((- kHarm(l+1,1)*dfhbm_dudot.*cl).*ck/Nfft,3);
                                    Judot((2*k-0)*NDof+(1:NDof),(2*l-1)*NInput+(1:NInput)) = -2*sum((- kHarm(l+1,1)*dfhbm_dudot.*sl).*sk/Nfft,3);
                                    Judot((2*k-0)*NDof+(1:NDof),(2*l-0)*NInput+(1:NInput)) = -2*sum((- kHarm(l+1,1)*dfhbm_dudot.*cl).*sk/Nfft,3);
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
                dfhbm_dxddot = hbm_derivatives('nl' ,'xddot',t,x,xdot,xddot,u,udot,uddot,xalg,f_hbm,hbm,problem,w0);
                dfalg_dxddot = hbm_derivatives('alg','xddot',t,x,xdot,xddot,u,udot,uddot,xalg,f_alg,hbm,problem,w0);
                
                if isfield(problem,'jacobXddot')
                    Jxddot = feval(problem.jacobXddot,dfhbm_dxddot,dfalg_dxddot,hbm,problem,w0);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Jxx = sum(hbm.nonlin.hbm.Jxddot{1}.*dfhbm_dxddot(ijacobx,ijacobx,:),3);
                            Jxa = zeros(NRetain,NAlg*prod(Nfft));
                            Jax = catmat(hbm.nonlin.alg.Jxddot{1}.*dfalg_dxddot(:,ijacobx,:),1);
                            Jaa = zeros(NAlg*prod(Nfft));
                            Jxddot = [Jxx Jxa;
                                        Jax Jaa];
                        case 'sum'
                            theta = 2*pi/Nfft*(0:(Nfft-1));
                            
                            Jxddot = zeros(NDof*NComp);
                            for l = 1:(NFreq-1)
                                [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                Jxddot(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum((-kHarm(l+1,1)^2*dfhbm_dxddot.*cl)/Nfft,3);
                                Jxddot(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum(( kHarm(l+1,1)^2*dfhbm_dxddot.*sl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                [ck, sk] = cosAndSin(kHarm(k+1,1)*theta);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                    Jxddot((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum((- kHarm(l+1,1)^2*dfhbm_dxddot.*cl).*ck/Nfft,3);
                                    Jxddot((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((  kHarm(l+1,1)^2*dfhbm_dxddot.*sl).*ck/Nfft,3);
                                    Jxddot((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum((- kHarm(l+1,1)^2*dfhbm_dxddot.*cl).*sk/Nfft,3);
                                    Jxddot((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((  kHarm(l+1,1)^2*dfhbm_dxddot.*sl).*sk/Nfft,3);
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
                dfhbm_duddot = hbm_derivatives('nl' ,'uddot',t,x,xdot,xddot,u,udot,uddot,xalg,f_hbm,hbm,problem,w0);
                dfalg_duddot = hbm_derivatives('alg','uddot',t,x,xdot,xddot,u,udot,uddot,xalg,f_alg,hbm,problem,w0);
                
                if isfield(problem,'jacobUddot')
                    Juddot = feval(problem.jacobUddot,dfhbm_duddot,dfalg_duddot,hbm,problem,w0);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Jxx = sum(hbm.nonlin.hbm.Juddot{1}.*dfhbm_duddot(ijacobx,ijacobu,:),3);
                            Jxa = zeros(NRetain,NAlg*prod(Nfft));
                            Jax = catmat(hbm.nonlin.alg.Juddot{1}.*dfalg_duddot(:,ijacobu,:),1);
                            Jaa = zeros(NAlg*prod(Nfft));
                            Juddot = [Jxx Jxa;
                                        Jax Jaa];
                        case 'sum'
                            theta = 2*pi/Nfft*(0:(Nfft-1));
                            
                            Juddot = zeros(NDof*NComp);
                            for l = 1:(NFreq-1)
                                [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                Juddot(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum((-kHarm(l+1,1)^2*dfhbm_duddot.*cl)/Nfft,3);
                                Juddot(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum(( kHarm(l+1,1)^2*dfhbm_duddot.*sl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                [ck, sk] = cosAndSin(kHarm(k+1,1)*theta);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                    Juddot((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum((- kHarm(l+1,1)^2*dfhbm_duddot.*cl).*ck/Nfft,3);
                                    Juddot((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((  kHarm(l+1,1)^2*dfhbm_duddot.*sl).*ck/Nfft,3);
                                    Juddot((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum((- kHarm(l+1,1)^2*dfhbm_duddot.*cl).*sk/Nfft,3);
                                    Juddot((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((  kHarm(l+1,1)^2*dfhbm_duddot.*sl).*sk/Nfft,3);
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
                dfhbm_dxdot = hbm_derivatives('nl' ,'xdot',t,x,xdot,xddot,u,udot,uddot,xalg,f_hbm,hbm,problem,w0);
                dfalg_dxdot = hbm_derivatives('alg','xdot',t,x,xdot,xddot,u,udot,uddot,xalg,f_alg,hbm,problem,w0);
                
                if isfield(problem,'floquet1xdot')
                    D1 = feval(problem.floquet1xdot,dfhbm_dxdot,dfalg_dxdot,hbm,problem,w0);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Jxx = sum(hbm.nonlin.hbm.Jx.*dfhbm_dxdot(ijacobx,ijacobx,:),3);
                            Jxa = zeros(NRetain,NAlg*prod(Nfft));
                            Jax = catmat(hbm.nonlin.alg.Jx.*dfalg_dxdot(:,ijacobx,:),1);
                            Jaa = zeros(NAlg*prod(Nfft));
                            D1 = [Jxx Jxa;
                                  Jax Jaa];
                        case 'sum'
                            theta = 2*pi/Nfft*(0:(Nfft-1));
                            
                            D1 = zeros(NDof*NComp);
                            D1(1:NDof,1:NDof) = mean(dfhbm_dxdot,3);
                            for l = 1:(NFreq-1)
                                [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                D1(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum(( dfhbm_dxdot.*cl)/Nfft,3);
                                D1(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum((-dfhbm_dxdot.*sl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                [ck, sk] = cosAndSin(kHarm(k+1,1)*theta);
                                D1((2*k-1)*NDof+(1:NDof),1:NDof) =  2*sum(dfhbm_dxdot.*ck/Nfft,3);
                                D1((2*k-0)*NDof+(1:NDof),1:NDof) = -2*sum(dfhbm_dxdot.*sk/Nfft,3);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                    D1((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum(( dfhbm_dxdot.*cl).*ck/Nfft,3);
                                    D1((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((-dfhbm_dxdot.*sl).*ck/Nfft,3);
                                    D1((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum(( dfhbm_dxdot.*cl).*sk/Nfft,3);
                                    D1((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((-dfhbm_dxdot.*sl).*sk/Nfft,3);
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
                dfhbm_dxddot = hbm_derivatives('nl' ,'xddot',t,x,xdot,xddot,u,udot,uddot,xalg,f_hbm,hbm,problem,w0);
                dfalg_dxddot = hbm_derivatives('alg','xddot',t,x,xdot,xddot,u,udot,uddot,xalg,f_alg,hbm,problem,w0);
                
                if isfield(problem,'floquet1xddot')
                    D1 = feval(problem.floquet1xddot,dfhbm_dxddot,dfalg_dxddot,hbm,problem,w0);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Jxx = sum(hbm.nonlin.hbm.Jxdot{1}.*dfhbm_dxddot(ijacobx,ijacobx,:),3);
                            Jxa = zeros(NRetain,NAlg*prod(Nfft));
                            Jax = catmat(hbm.nonlin.alg.Jxdot{1}.*dfalg_dxddot(:,ijacobx,:),1);
                            Jaa = zeros(NAlg*prod(Nfft));
                            D1 = [Jxx Jxa;
                                  Jax Jaa];
                        case 'sum'
                            theta = 2*pi/Nfft*(0:(Nfft-1));
                            
                            D1 = zeros(NDof*NComp);
                            for l = 1:(NFreq-1)
                                [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                D1(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum((-kHarm(l+1,1)*dfhbm_dxddot.*sl)/Nfft,3);
                                D1(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum((-kHarm(l+1,1)*dfhbm_dxddot.*cl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                [ck, sk] = cosAndSin(kHarm(k+1,1)*theta);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                    D1((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum((-kHarm(l+1,1)*dfhbm_dxddot.*sl).*ck/Nfft,3);
                                    D1((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((-kHarm(l+1,1)*dfhbm_dxddot.*cl).*ck/Nfft,3);
                                    D1((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum((-kHarm(l+1,1)*dfhbm_dxddot.*sl).*sk/Nfft,3);
                                    D1((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((-kHarm(l+1,1)*dfhbm_dxddot.*cl).*sk/Nfft,3);
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
                dfhbm_dxddot = hbm_derivatives('nl' ,'xddot',t,x,xdot,xddot,u,udot,uddot,xalg,f_hbm,hbm,problem,w0);
                dfalg_dxddot = hbm_derivatives('alg','xddot',t,x,xdot,xddot,u,udot,uddot,xalg,f_alg,hbm,problem,w0);
                
                if isfield(problem,'floquet2')
                    D2 = feval(problem.floquet2,dfhbm_dxddot,dfalg_dxddot,hbm,problem,w0);
                else
                    switch hbm.options.jacob_method
                        case 'mat'
                            Jxx = sum(hbm.nonlin.hbm.Jx.*dfhbm_dxddot(ijacobx,ijacobx,:),3);
                            Jxa = zeros(NRetain,NAlg*prod(Nfft));
                            Jax = catmat(hbm.nonlin.alg.Jx.*dfalg_dxddot(:,ijacobx,:),1);
                            Jaa = zeros(NAlg*prod(Nfft));
                            D2 = [Jxx Jxa;
                                  Jax Jaa];
                        case 'sum'
                            theta = 2*pi/Nfft*(0:(Nfft-1));
                            
                            D2 = zeros(NDof*NComp);
                            D2(1:NDof,1:NDof) = mean(dfhbm_dxddot,3);
                            for l = 1:(NFreq-1)
                                [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                D2(1:NDof,(2*l-1)*NDof+(1:NDof)) =  sum(( dfhbm_dxddot.*cl)/Nfft,3);
                                D2(1:NDof,(2*l-0)*NDof+(1:NDof)) =  sum((-dfhbm_dxddot.*sl)/Nfft,3);
                            end
                            for k = 1:(NFreq-1)
                                [ck, sk] = cosAndSin(kHarm(k+1,1)*theta);
                                D2((2*k-1)*NDof+(1:NDof),1:NDof) =  2*sum(dfhbm_dxddot.*ck/Nfft,3);
                                D2((2*k-0)*NDof+(1:NDof),1:NDof) = -2*sum(dfhbm_dxddot.*sk/Nfft,3);
                                for l = 1:(NFreq-1)
                                    [cl, sl] = cosAndSin(kHarm(l+1,1)*theta);
                                    D2((2*k-1)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) =  2*sum(( dfhbm_dxddot.*cl).*ck/Nfft,3);
                                    D2((2*k-1)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) =  2*sum((-dfhbm_dxddot.*sl).*ck/Nfft,3);
                                    D2((2*k-0)*NDof+(1:NDof),(2*l-1)*NDof+(1:NDof)) = -2*sum(( dfhbm_dxddot.*cl).*sk/Nfft,3);
                                    D2((2*k-0)*NDof+(1:NDof),(2*l-0)*NDof+(1:NDof)) = -2*sum((-dfhbm_dxddot.*sl).*sk/Nfft,3);
                                end
                            end
                    end
                    D2 = D2(iRetain,iRetain);
                end
            end
            varargout{o} = D2;
           
    end
end

function [cl, sl] = cosAndSin(ph)
cl = cos(ph); sl = sin(ph);
