function problem = test_params
rng(5);
P.k1 = 10+rand(1)*0;
P.k2 = 10+rand(1)*0;
P.k3 = 5+rand(1)*0;
P.m1 = 1;
P.m2 = 2;
P.c1 = 0.05;
P.c2 = 0.05;
P.c3 = 0.05;

P.cnl = 0*[0 0 0.005]';
P.knl = [0 0 0.004]';
P.n = 3;

P.R = [1 0;
      -1 1;
       0 1];

%forcing
P.f0 = 0;
P.f = 2*rand(1)*exp(1i*pi*rand(1));
P.f2 = 0.5*rand(1)*exp(1i*pi*rand(1));

P.iDof = 1;
P.iInput = 1;

problem.K = P.k1*[1 0;
                  0 0];
problem.K = problem.K + P.k2*[1 -1;
             -1  1];    
problem.K = problem.K + P.k3*[0 0;
             0  1];    
problem.C = P.c1*[1 0;
          0 0];
problem.C = problem.C + P.c2*[1 -1;
             -1  1];
problem.C = problem.C + P.c3*[0 0;
             0  1];    
problem.M = diag([P.m1 P.m2]);
         
problem.F0 = zeros(2,1);

problem.Mu = zeros(2,1);
problem.Cu = zeros(2,1);
problem.Ku = [0;1];

problem.model = @test_model;
problem.excite = @test_excite;

problem.P = P;