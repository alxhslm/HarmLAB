function varargout = get_time_series(hbm,w0,X)
NHarm = hbm.harm.NHarm(1);
Nfft = hbm.harm.Nfft(1);

%unpack the inputs
w = (0:NHarm)'*w0;
t = (0:Nfft-1)'/Nfft*2*pi/w0;

%compute the fourier coefficients of the derivatives
Wx = repmat(1i*w,1,size(X,2));

switch hbm.options.aft_method
    case 'fft'
        %create the time series from the fourier series
        x = freq2time(X,NHarm,Nfft);
        if nargout > 2
            Xdot  = X.*Wx;
            xdot  = freq2time(Xdot,NHarm,Nfft);
            if nargout > 3
                Xddot = Xdot.*Wx;
                xddot = freq2time(Xddot,NHarm,Nfft);
            end
        end
        
    case 'mat'
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
end

varargout{1} = t;
varargout{2} = x;
if nargout >2
    varargout{3} = xdot;
    if nargout > 3
        varargout{4} = xddot;
    end
end