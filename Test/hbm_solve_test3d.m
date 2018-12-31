function hbm_solve_test3d
problem = test_params;
NDof = problem.NDof;
NAlg = problem.NAlg;

hbm.disp_dependence = true;
hbm.vel_dependence = true;
hbm.freq_dependence = false;

hbm.harm.rFreqRatio = 1;
hbm.harm.NHarm = 5;
hbm.harm.Nfft  = 16;

if 1
    hbm.harm.rFreqRatio(end+1) = 1.376;
    hbm.harm.NHarm(end+1) = 3;
    hbm.harm.Nfft(end+1) = 16;
end

hbm.options.disp_dependence = true;
hbm.options.vel_dependence  = true;
hbm.options.freq_dependence = false;

hbm.options.cont_method = 'coco';
hbm.options.aft_method = 'mat';
hbm.options.jacob_method = 'mat';

hbm = setuphbm(hbm,problem);

% X = rand(5,1);
% x = real(hbm.nonlin.IFFT*X);
% X2 = hbm.nonlin.FFT*x;
% x = reshape(x,hbm.harm.Nfft(1),hbm.harm.Nfft(2));
% surf(x)
% hbm_opt = hbm;

% d =hbm2.aft.FFT  - hbm.nonlin.FFT;
% max(abs(d(:)))
% d =hbm2.aft.IFFT - hbm.nonlin.IFFT;
% max(abs(d(:)))
% d =hbm2.nonlin.Jx  - hbm.nonlin.Jx;
% max(abs(d(:)))
% d=hbm2.nonlin.Jxdot{1}  - hbm.nonlin.Jxdot;
% max(abs(d(:)))
% d=hbm2.nonlin.Ju - hbm.nonlin.Ju;
% max(abs(d(:)))
% d=hbm2.nonlin.Judot{1}  - hbm.nonlin.Judot;
% max(abs(d(:)))


w0 = 2.05;
A = 0.5;

hbm.options.aft_method = 'fft';
hbm.options.jacob_method = 'mat';
hbm = setuphbm(hbm,problem);
tic;
[X1,XAlg1,U,F] = hbm_solve(hbm,problem,w0,A,[]);
toc;
[t1,x1,xalg1,xdot1] = get_time_series3d(hbm,w0,X1,XAlg1);

hbm.options.aft_method = 'fft';
hbm.options.jacob_method = 'sum';
hbm = setuphbm(hbm,problem);
tic;
x0 = [packdof(X1);reshape(XAlg1,1,[])'];
[X2,XAlg2,U,F] = hbm_solve(hbm,problem,w0,A,x0);
toc;
[t2,x2,xalg2,xdot2] = get_time_series3d(hbm,w0,X2,XAlg2);

hbm.options.aft_method = 'mat';
hbm.options.jacob_method = 'mat';
hbm = setuphbm(hbm,problem);
tic;
x0 = [packdof(X1);reshape(XAlg2,1,[])'];
[X3,XAlg3,U,F] = hbm_solve(hbm,problem,w0,A,x0);
toc;
[t3,x3,xalg3,xdot3] = get_time_series3d(hbm,w0,X3,XAlg3);

hbm.options.aft_method = 'mat';
hbm.options.jacob_method = 'sum';
hbm = setuphbm(hbm,problem);
tic;
x0 = [packdof(X1);reshape(XAlg3,1,[])'];
[X4,XAlg4,U,F] = hbm_solve(hbm,problem,w0,A,x0);
toc;
[t4,x4,xalg4,xdot4] = get_time_series3d(hbm,w0,X4,XAlg4);

% tic;
y0 = [x4(1,:) xdot4(1,:) xalg4(1,:)];
fun = @(t,y)test_odefun(t,y,w0*hbm.harm.rFreqRatio,U,hbm,problem);
M = blkdiag(eye(NDof),problem.M,zeros(NAlg));
options = odeset('Mass',M,'Vectorized','on','MassSingular','yes','RelTol',1E-12,'AbsTol',1E-12);
[t5,y] = ode15s(fun,[t4(1) t4(end)],y0,options);
tRun(5) = toc;

x5 = y(:,1:NDof);
xdot5 = y(:,NDof+(1:NDof));
xalg5 = y(:,2*NDof+(1:NAlg));

% X5 = time2freq3d(x5,hbm.harm.NHarm,hbm.iSub,hbm.harm.Nfft);

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

subplot(3,2,3)
plot(t1,xalg1','-',t2,xalg2','-',t3,xalg3','-',t4,xalg4','-',t5,xalg5','o-');
ylabel('x_{alg}');

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