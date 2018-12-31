function G = test_obj(X,U,F,hbm,problem,w0)
P = problem.P;
H = X(2,P.iDof)./U(2,P.iInput);
% G = -(angle(H) + pi/2).^2;
G = abs(H);