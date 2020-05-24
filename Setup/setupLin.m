function lin = setupLin(harm,problem)
%% state-jacobian
[lin.Ak,lin.Ac,lin.Am,lin.Ax] = linear_jacobian(harm,problem.K,problem.C,problem.M,problem.G);

%% input-jacobian
[lin.Bk,lin.Bc,lin.Bm,lin.Bx] = linear_jacobian(harm,problem.Ku,problem.Cu,problem.Mu);

%% constant terms
lin.b = [problem.F0;
        zeros(problem.NDof*(harm.NComp-1),1)];

%% floquet multipliers
[floquet.D1xdot,floquet.D1xddot] = linear_jacobian(harm,problem.C,problem.M,0*problem.M);
floquet.D2 = kron(eye(harm.NComp),problem.M);

lin.floquet = floquet;