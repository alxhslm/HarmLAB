function nonlin = setupNonlin(harm,problem)
Nfft = harm.Nfft;

%% fft matrices
theta1 = 2*pi/Nfft(1)*(0:(Nfft(1)-1));
theta2 = 2*pi/Nfft(2)*(0:(Nfft(2)-1));
[theta1,theta2] = ndgrid(theta1,theta2);
theta0 = [theta1(:)'; theta2(:)'];

nonlin.FFT = zeros(harm.NFreq,Nfft(1)*Nfft(2));
nonlin.FFT(1,:) = 1/(Nfft(1)*Nfft(2));
for k = 2:harm.NFreq
    nonlin.FFT(k,:) =  2/(Nfft(1)*Nfft(2))*exp(-1i*harm.kHarm(k,:)*theta0);
end

nonlin.IFFT = zeros(Nfft(1)*Nfft(2),harm.NFreq);
nonlin.IFFT(:,1) = 1;
for k = 2:harm.NFreq
    nonlin.IFFT(:,k) = exp(1i*harm.kHarm(k,:)*theta0);
end

%% non-linear derivative matrices 
%work out the constituent jacobians first
aft = get_aft_jacobians(harm);

%now make the specifis matrices we will need
nonlin.hbm = get_hbm_matrices(problem,harm,aft);

function aft = get_aft_jacobians(harm)
aft.fft.J = aft_jacobian_fft(harm);
[aft.ifft.J,aft.ifft.Jdot,aft.ifft.Jddot] = aft_jacobian_ifft(harm);

%total jacobian
aft.J = mtimesx(aft.fft.J,aft.ifft.J);
for i = 1:2
    aft.Jdot{i}  = mtimesx(aft.fft.J,harm.rFreqBase(i)*aft.ifft.Jdot{i});
    aft.Jddot{i} = mtimesx(aft.fft.J,harm.rFreqBase(i)^2*aft.ifft.Jddot{i});
end
aft.Jddot{3} = mtimesx(aft.fft.J,prod(harm.rFreqBase)*aft.ifft.Jddot{3});

function hbm = get_hbm_matrices(problem,harm,aft)

%hbm
hbm.Jfft = resize_jacobian(repmat(aft.fft.J,1,size(aft.fft.J,1)),problem.NDof,problem.NNL); 

hbm.Jifft = resize_jacobian(repmat(aft.ifft.J,size(aft.ifft.J,2),1),problem.NDof,problem.NNL); 

for i = 1:2
    hbm.Jdotifft{i} = resize_jacobian(repmat(harm.rFreqBase(i)*aft.ifft.Jdot{i},size(aft.ifft.J,2),1),problem.NDof,problem.NNL); 
end

%% States
hbm.Jx     = resize_jacobian(aft.J,problem.NDof,problem.NNL); 
hbm.Jxdot  = resize_jacobian(aft.Jdot,problem.NDof,problem.NNL); 
hbm.Jxddot  = resize_jacobian(aft.Jddot,problem.NDof,problem.NNL);

%% Inputs
hbm.Ju     = resize_jacobian(aft.J,problem.NDof,problem.NInput); 
hbm.Judot  = resize_jacobian(aft.Jdot,problem.NDof,problem.NInput); 
hbm.Juddot  = resize_jacobian(aft.Jddot,problem.NDof,problem.NInput); 


%indices for creating the jacobian
hbm.ijacobx = repmat((1:problem.NDof)',harm.NComp,1);
hbm.ijacobu = repmat((1:problem.NInput)',harm.NComp,1);
hbm.ijacobxnl = repmat(problem.iNL(:),harm.NComp,1);

function Ju = resize_jacobian(ju,NOutput,NInput)
if ~iscell(ju)
    ju = {ju};
end
for k = 1:length(ju)
    Ju{k} = zeros(NOutput*size(ju{k},1),NInput*size(ju{k},2),size(ju{k},3));
    for i = 1:size(ju{k},3)
        Ju{k}(:,:,i) = kron(ju{k}(:,:,i),ones(NOutput,NInput));
    end
end
if length(Ju) == 1
    Ju = Ju{1};
end