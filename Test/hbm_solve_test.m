function hbm_solve_test
problem = test_params;
NDof = problem.NDof;
NAlg = problem.NAlg;

hbm.harm.NHarm = 2;
hbm.harm.Nfft = 32;

hbm.options.bUseStandardHBM = true;
hbm.options.disp_dependence = true;
hbm.options.vel_dependence  = true;
hbm.options.freq_dependence = false;

hbm.options.bAnalyticalDerivs = true;
hbm.options.cont_method = 'coco';
hbm.options.aft_method = 'mat';
hbm.options.jacob_method = 'mat';

omega = sqrt(eig(problem.K,problem.M));
w0 = 4;
A = 100;

hbm = setuphbm(hbm,problem);

hbm.options.aft_method = 'mat';
hbm.options.jacob_method = 'mat';
hbm = setuphbm(hbm,problem);
tic;
[X1,xAlg1,U,F] = hbm_solve(hbm,problem,w0,A,[]);
[t1,x1,xdot1] = get_time_series(hbm,w0,X1);
tRun(1) = toc;

hbm.options.aft_method = 'mat';
hbm.options.jacob_method = 'sum';
hbm = setuphbm(hbm,problem);
tic;
x0 = [packdof(X1);reshape(xAlg1,1,[])'];
[X2,xAlg2,U,F] = hbm_solve(hbm,problem,w0,A,x0);
[t2,x2,xdot2] = get_time_series(hbm,w0,X2);
tRun(2) = toc;

hbm.options.aft_method = 'fft';
hbm.options.jacob_method = 'mat';
hbm = setuphbm(hbm,problem);
tic;
x0 = [packdof(X2);reshape(xAlg2,1,[])'];
[X3,xAlg3,U,F] = hbm_solve(hbm,problem,w0,A,x0);
[t3,x3,xdot3] = get_time_series(hbm,w0,X3);
tRun(3) = toc;

hbm.options.aft_method = 'fft';
hbm.options.jacob_method = 'sum';
hbm = setuphbm(hbm,problem);
tic;
x0 = [packdof(X3);reshape(xAlg3,1,[])'];
[X4,xAlg4,U,F] = hbm_solve(hbm,problem,w0,A,x0);
[t4,x4,xdot4] = get_time_series(hbm,w0,X4);
tRun(4) = toc;

tic;
y0 = [x4(1,:) xdot4(1,:) xAlg4(1,:)];
fun = @(t,y)test_odefun(t,y,w0,U,hbm,problem);
M = blkdiag(eye(NDof),problem.M,zeros(NAlg));
options = odeset('Mass',M,'Vectorized','on','MassSingular','yes','RelTol',1E-12,'AbsTol',1E-12);
[t5,y] = ode15s(fun,[t4(1) t4(end)],y0,options);
tRun(5) = toc;

x5 = y(:,1:NDof);
xdot5 = y(:,NDof+(1:NDof));
xAlg5 = y(:,2*NDof+(1:NAlg));
% for i = 1:size(x2,1)
%     xInt(i,1) = fsolve(@(z)getnthoutput(2,2,@force,P,x2(i,:)',z),0,optimoptions('fsolve','Display','off'));
% end
% x2 = [x2 xInt];
% xdot2 = [xdot2, adiff(t2,xInt,1)];
% toc;

figure
subplot(3,2,1)
plot(t1,x1(:,1:NDof)','-',t2,x2(:,1:NDof)','-',t3,x3(:,1:NDof)','-',t4,x4(:,1:NDof)','-',t5,x5(:,1:NDof)','o-');
ylabel('x');

if problem.NAlg>0
    subplot(3,2,3)
    plot(t1,xAlg1','-',t2,xAlg2','-',t3,xAlg3','-',t4,xAlg4','-',t5,xAlg5','o-');
    ylabel('x_{alg}');
end

subplot(3,2,5)
plot(t1,xdot1(:,1:NDof)','-',t2,xdot2(:,1:NDof)','-',t3,xdot3(:,1:NDof)','-',t4,xdot4(:,1:NDof)','-',t5,xdot5(:,1:NDof)','o-');
ylabel('xdot');
xlabel('Time (s)');

names = {'HBM-mat/mat','HBM-mat/sum','HBM-fft/mat','HBM-fft/sum','ODE'};

subplot(1,2,2)
plot(x1(:,1:NDof),xdot1(:,1:NDof));
hold on
plot(x2(:,1:NDof),xdot2(:,1:NDof));
plot(x3(:,1:NDof),xdot3(:,1:NDof));
plot(x4(:,1:NDof),xdot4(:,1:NDof));
plot(x5(:,1:NDof),xdot5(:,1:NDof),'o-');
xlabel('x');
ylabel('xdot');
leg = {};
for i = 1:5
    leg = [leg cellsprintf([names{i} '(x%d)'],num2cell(1:NDof))];
end
legend(leg)
    
for i = 1:5
     fprintf('%s: %f s\n',names{i},tRun(i));
end