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

%% input-jacobian
[lin.Bk,lin.Bc,lin.Bm,lin.Bx] = linear_jacobian(harm,problem.Ku,problem.Cu,problem.Mu);

%retain only terms we need
for i = 1:length(lin.Bk)
    lin.Bk{i} = lin.Bk{i}(iRetain,:);
    lin.Bc{i} = lin.Bc{i}(iRetain,:);
    lin.Bm{i} = lin.Bm{i}(iRetain,:);
end
lin.Bx = lin.Bx(iRetain,:);

%% constant terms
lin.b = [problem.F0;
        zeros(problem.NDof*(harm.NComp-1),1)];
lin.b = lin.b(iRetain);

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