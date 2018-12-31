function varargout = get_time_series3d(hbm,w0,X,xalg,tspan)
if hbm.options.bUseStandardHBM
    [varargout{1:nargout}] = get_time_series(hbm,problem,w0,u,X);
    return;
end

NHarm = hbm.harm.NHarm;
Nfft = hbm.harm.Nfft;
kHarm  = hbm.harm.kHarm;
w0 = w0 * hbm.harm.rFreqRatio;

%unpack the inputs
w = kHarm(:,1)*w0(1) + kHarm(:,2)*w0(2);
t1 = (0:Nfft(1)-1)/Nfft(1)*2*pi/w0(1);
t2 = (0:Nfft(2)-1)/Nfft(2)*2*pi/w0(2);
[t1,t2] = ndgrid(t1,t2);
%compute the fourier coefficients of the derivatives
Wx = repmat(1i*w,1,size(X,2));

switch hbm.options.aft_method
    case 'fft'
        %create the time series from the fourier series
        x = freq2time3d(X,NHarm,hbm.harm.iSub,Nfft);
        if nargout > 2
            Xdot  = X.*Wx;
            xdot  = freq2time3d(Xdot,NHarm,hbm.harm.iSub,Nfft);
            if nargout > 3
                Xddot = Xdot.*Wx;
                xddot = freq2time3d(Xddot,NHarm,hbm.harm.iSub,Nfft);
            end
        end
        
    case 'mat'
        %create the time series from the fourier series
        x = real(hbm.nonlin.IFFT*X);
        if nargout > 2
            Xdot  = X.*Wx;
            xdot  = real(hbm.nonlin.IFFT*Xdot);
            if nargout > 3
                Xddot  = Xdot.*Wx;
                xddot = real(hbm.nonlin.IFFT*Xddot);
            end
        end
end

x = reshape(x,size(t1,1),size(t1,2),[]);
if nargout > 2
    xdot = reshape(xdot,size(t1,1),size(t1,2),[]);
    if nargout > 3
        xddot = reshape(xddot,size(t1,1),size(t1,2),[]);
    end
end

if nargin < 5 || isempty(tspan) %empty
    ti = linspace(0,2*pi/min(w0),max(Nfft));
elseif length(tspan) < 2 %duration
    ti = linspace(0,tspan,max(Nfft));
else
    ti = tspan;
end
    
if ti(end) > 2*pi/w0(1)% && hbm.harm.NHarm(1)>0
    N = ceil(ti(end)/t1(end,1));
    t0 = t1;
    for i = 2:N
        t1 = [t1; t0+2*pi*(i-1)/w0(1)];
    end
    x = repmat(x,N,1);
    t2 = repmat(t2,N,1);
    if nargout > 2
        xdot = repmat(xdot,N,1);
        if nargout > 3
            xddot = repmat(xddot,N,1);
        end
    end
    if nargin > 3 && ~isempty(xalg)
        xalg = repmat(xalg,N,1);
    end
end
if ti(end) > 2*pi/w0(2)% && hbm.harm.NHarm(2)>0
    N = ceil(ti(end)/t2(1,end));
    t0 = t2;
    for i = 2:N
        t2 = [t2  t0+2*pi*(i-1)/w0(2)];
    end
    x = repmat(x,1,N);
    t1 = repmat(t1,1,N);
    if nargout > 2
        xdot = repmat(xdot,1,N);
        if nargout > 3
            xddot = repmat(xddot,1,N);
        end
    end
    if nargin > 3 && ~isempty(xalg)
        xalg = repmat(xalg,1,N);
    end
end

if Nfft(2) > 1
    for i = 1:size(x,3)
        xi(:,i) = interp2(t2,t1,x(:,:,i),ti,ti);
        if nargout > 2
            xdoti(:,i) = interp2(t2,t1,xdot(:,:,i),ti,ti);
            if nargout > 3
                xddoti(:,i) = interp2(t2,t1,xddot(:,:,i),ti,ti);
            end
        end
    end
    if nargin > 3 && ~isempty(xalg)
        for i = 1:size(xalg,3)
            xalgi(:,i) = interp2(t2,t1,xalg(:,:,i),ti,ti);
        end
    end
else
    for i = 1:size(x,3)
        xi(:,i) = interp1(t1,x(:,:,i),ti);
        if nargout > 2
            xdoti(:,i) = interp1(t1,xdot(:,:,i),ti);
            if nargout > 3
                xddoti(:,i) = interp1(t1,xddot(:,:,i),ti);
            end
        end
        if nargin > 3 && ~isempty(xalg)
            for i = 1:size(xalg,3)
                xalgi(:,i) = interp1(t1,xalg(:,:,i),ti);
            end
        end
    end
end

varargout{1} = ti;
varargout{2} = xi;
if nargin > 3 && ~isempty(xalg)
    varargout{3} = xalgi;
end

if nargout >2
    varargout{4} = xdoti;
    if nargout > 3
        varargout{5} = xddoti;
    end
end