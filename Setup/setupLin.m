function lin = setupLin(harm,problem)
iRetain = harm.iRetain;

%% state-jacobian
[lin.Ak,lin.Ac,lin.Am,lin.Ax] = linear_jacobian(harm,problem.K,problem.C,problem.M);

%retain only terms we need
for i = 1:length(lin.Ak)
    lin.Ak{i} = lin.Ak{i}(iRetain,iRetain);
    lin.Ac{i} = lin.Ac{i}(iRetain,iRetain);
    lin.Am{i} = lin.Am{i}(iRetain,iRetain);
end
lin.Ax = lin.Ax(iRetain,iRetain);

%pad with zeros for alg constraints
for i = 1:length(lin.Ak)
    lin.Ak{i} = blkdiag(lin.Ak{i}, zeros(problem.NAlg*prod(harm.Nfft)));
    lin.Ac{i} = blkdiag(lin.Ac{i}, zeros(problem.NAlg*prod(harm.Nfft)));
    lin.Am{i} = blkdiag(lin.Am{i}, zeros(problem.NAlg*prod(harm.Nfft)));
end
lin.Ax = blkdiag(lin.Ax, zeros(problem.NAlg*prod(harm.Nfft)));

%% input-jacobian
[lin.Bk,lin.Bc,lin.Bm,lin.Bx] = linear_jacobian(harm,problem.Ku,problem.Cu,problem.Mu);

%retain only terms we need
for i = 1:length(lin.Bk)
    lin.Bk{i} = lin.Bk{i}(iRetain,:);
    lin.Bc{i} = lin.Bc{i}(iRetain,:);
    lin.Bm{i} = lin.Bm{i}(iRetain,:);
end
lin.Bx = lin.Bx(iRetain,:);

for i = 1:length(lin.Bk)
    lin.Bk{i} = [lin.Bk{i}; zeros(problem.NAlg*prod(harm.Nfft),harm.NComp*problem.NInput)];
    lin.Bc{i} = [lin.Bc{i}; zeros(problem.NAlg*prod(harm.Nfft),harm.NComp*problem.NInput)];
    lin.Bm{i} = [lin.Bm{i}; zeros(problem.NAlg*prod(harm.Nfft),harm.NComp*problem.NInput)];
end
lin.Bx = [lin.Bx; zeros(problem.NAlg*prod(harm.Nfft),size(lin.Bx,2))];

%% constant terms
lin.b = [problem.F0;
        zeros(problem.NDof*(harm.NComp-1),1)];
lin.b = [lin.b(iRetain);
    zeros(problem.NAlg*prod(harm.Nfft),1)];

%% floquet multipliers
[floquet.D1xdot,floquet.D1xddot] = linear_jacobian(harm,problem.C,problem.M,0*problem.M);
floquet.D1xdot = floquet.D1xdot{1};
floquet.D2 = kron(eye(harm.NComp),problem.M);

%trim to size
for i = 1:2
    floquet.D1xddot{i} = floquet.D1xddot{i}(iRetain,iRetain);
end
floquet.D1xdot = floquet.D1xdot(iRetain,iRetain);
floquet.D2 = floquet.D2(iRetain,iRetain);

lin.floquet = floquet;