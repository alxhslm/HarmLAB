function [Ju,Judot] = aft_jacobian(harm,NOutput,NInput)
Nfft = harm.Nfft;
NComp = harm.NComp;

theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
[theta1,theta2] = ndgrid(theta1,theta2);
theta0 = [theta1(:)'; theta2(:)'];

theta1 = permute(theta0(1,:),[1 3 2]);
theta2 = permute(theta0(2,:),[1 3 2]);

ju = zeros(NComp,NComp,Nfft(1)*Nfft(2));
ju(1,1,:) = 1/(Nfft(1)*Nfft(2));
for l = 1:(harm.NFreq-1)
    ph = harm.kHarm(l+1,1)*theta1 + harm.kHarm(l+1,2)*theta2;
    cl = cos(ph); sl = sin(ph);
    ju(1,2*l,:) =  cl/(Nfft(1)*Nfft(2));
    ju(1,2*l+1,:) = -sl/(Nfft(1)*Nfft(2));
end
for k = 1:(harm.NFreq-1)
    ph = harm.kHarm(k+1,1)*theta1 + harm.kHarm(k+1,2)*theta2;
    ck = cos(ph); sk = sin(ph);
    ju(2*k,1,:) =  2*ck/(Nfft(1)*Nfft(2));
    ju(2*k+1,1,:) = -2*sk/(Nfft(1)*Nfft(2));
    for l = 1:(harm.NFreq-1)
        ph = harm.kHarm(l+1,1)*theta1 + harm.kHarm(l+1,2)*theta2;
        cl = cos(ph); sl = sin(ph);
        ju(2*k,2*l,:) =  2*cl.*ck/(Nfft(1)*Nfft(2));
        ju(2*k,2*l+1,:) = -2*sl.*ck/(Nfft(1)*Nfft(2));
        ju(2*k+1,2*l,:) = -2*cl.*sk/(Nfft(1)*Nfft(2));
        ju(2*k+1,2*l+1,:) =  2*sl.*sk/(Nfft(1)*Nfft(2));
    end
end
Ju = resize_jacobian(ju,NInput,NOutput);

Judot = cell(1,2);
for n = 1:2
    judot = zeros(NComp,NComp,Nfft(1)*Nfft(2));
    for l = 1:(harm.NFreq-1)
        ph = harm.kHarm(l+1,1)*theta1 + harm.kHarm(l+1,2)*theta2;
        cl = cos(ph); sl = sin(ph);
        judot(1,2*l,:) = -harm.kHarm(l+1,n)*sl/(Nfft(1)*Nfft(2));
        judot(1,2*l+1,:) = -harm.kHarm(l+1,n)*cl/(Nfft(1)*Nfft(2));
    end
    for k = 1:(harm.NFreq-1)
        ph = harm.kHarm(k+1,1)*theta1 + harm.kHarm(k+1,2)*theta2;
        ck = cos(ph); sk = sin(ph);
        judot(2*k,1,:) =  0;
        judot(2*k+1,1,:) =  0;
        for l = 1:(harm.NFreq-1)
            ph = harm.kHarm(l+1,1)*theta1 + harm.kHarm(l+1,2)*theta2;
            cl = cos(ph); sl = sin(ph);
            judot(2*k,2*l,:) = -2*harm.kHarm(l+1,n).*sl.*ck/(Nfft(1)*Nfft(2));
            judot(2*k,2*l+1,:) = -2*harm.kHarm(l+1,n).*cl.*ck/(Nfft(1)*Nfft(2));
            judot(2*k+1,2*l,:) =  2*harm.kHarm(l+1,n).*sl.*sk/(Nfft(1)*Nfft(2));
            judot(2*k+1,2*l+1,:) =  2*harm.kHarm(l+1,n).*cl.*sk/(Nfft(1)*Nfft(2));
        end
    end
    Judot{n} = resize_jacobian(judot,NInput,NOutput);
end

function Ju = resize_jacobian(ju,NInput,NOutput)
Ju = zeros(NOutput*size(ju,1),NInput*size(ju,2),size(ju,3));
for i = 1:size(ju,3)
    Ju(:,:,i) = kron(ju(:,:,i),ones(NOutput,NInput));
end