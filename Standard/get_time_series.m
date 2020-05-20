function varargout = get_time_series(hbm,w0,X,tspan)
NHarm = hbm.harm.NHarm(1);
Nfft = hbm.harm.Nfft(1);

%unpack the inputs
w = (0:NHarm)'*w0;
t = (0:Nfft-1)'/Nfft*2*pi/w0;

%compute the fourier coefficients of the derivatives
Wx = repmat(1i*w,1,size(X,2));

%create the time series from the fourier series
x  = real(hbm.nonlin.IFFT*X);
if nargout > 2
    Xdot  = X.*Wx;
    xdot  = real(hbm.nonlin.IFFT*Xdot);
    if nargout > 3
        Xddot  = Xdot.*Wx;
        xddot = real(hbm.nonlin.IFFT*Xddot);
    end
end

if nargin < 5 || isempty(tspan) %empty
    ti = linspace(0,2*pi/min(w0),Nfft);
elseif length(tspan) < 2 %duration
    ti = linspace(0,tspan,Nfft);
else
    ti = tspan;
end

for i = 1:size(x,2)
    xi(:,i) = interp1(t,x(:,i),ti);
    if nargout > 2
        xdoti(:,i) = interp1(t,xdot(:,i),ti);
        if nargout > 3
            xddoti(:,i) = interp1(t,xddot(:,i),ti);
        end
    end
end

varargout{1} = ti;
varargout{2} = xi;
if nargout >2
    varargout{3} = xdoti;
    if nargout > 3
        varargout{4} = xddoti;
    end
end