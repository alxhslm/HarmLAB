function [Ju,Judot,Juddot] = aft_jacobian_ifft(harm)
Nfft = harm.Nfft;
NComp = harm.NComp;

theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
[theta1,theta2] = ndgrid(theta1,theta2);
theta0 = [theta1(:)'; theta2(:)'];

Ju = zeros(1,NComp,Nfft(1)*Nfft(2));
Ju(1,1,:) = 1;
for l = 1:(harm.NFreq-1)
    ph = permute(harm.kHarm(l+1,:)*theta0,[1 3 2]);
    cl = cos(ph); sl = sin(ph);
    Ju(1,2*l,:)   =  cl;
    Ju(1,2*l+1,:) = -sl;
end

Judot = cell(1,2);
for n = 1:2
    Judot{n} = zeros(1,NComp,Nfft(1)*Nfft(2));
    for l = 1:(harm.NFreq-1)
        ph = permute(harm.kHarm(l+1,:)*theta0,[1 3 2]);
        cl = cos(ph); sl = sin(ph);
        Judot{n}(1,2*l,:)   = -harm.kHarm(l+1,n)*sl;
        Judot{n}(1,2*l+1,:) = -harm.kHarm(l+1,n)*cl;
    end
end

Juddot = cell(1,3);
for n = 1:2
    Juddot{n} = zeros(1,NComp,Nfft(1)*Nfft(2));
    for l = 1:(harm.NFreq-1)
        ph = permute(harm.kHarm(l+1,:)*theta0,[1 3 2]);
        cl = cos(ph); sl = sin(ph);
        Juddot{n}(1,2*l,:)   = -harm.kHarm(l+1,n)^2*cl;
        Juddot{n}(1,2*l+1,:) =  harm.kHarm(l+1,n)^2*sl;
    end
end

n = 3;
Juddot{n} = zeros(1,NComp,Nfft(1)*Nfft(2));
for l = 1:(harm.NFreq-1)
    ph = permute(harm.kHarm(l+1,:)*theta0,[1 3 2]);
    cl = cos(ph); sl = sin(ph);
    Juddot{n}(1,2*l,:)   = -prod(harm.kHarm(l+1,:))*cl;
    Juddot{n}(1,2*l+1,:) =  prod(harm.kHarm(l+1,:))*sl;
end