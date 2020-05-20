function lambda = hbm_floquet(hbm,problem,w,U,X)
if ~isvector(X)
    x = packdof(X);
    u = packdof(U);
else
    x = X;
    u = U;
end

if any(isnan(x) | isinf(x))
    lambda = NaN(hbm.harm.NComp*problem.NDof,1);
    return
end

[A,B] = floquetMatrices(hbm,problem,w,u,x);
lambda = eig(A,B,'vector');
[~,iSort] = sort(abs(imag(lambda)));
lambda = lambda(iSort);


function [A,B] = floquetMatrices(hbm,problem,w,u,x)
hbm.bIncludeNL = 1;
NPts = size(u,2);

A = zeros(2*size(x,1),2*size(x,1),NPts);
B = zeros(2*size(x,1),2*size(x,1),NPts);
for i = 1:NPts
    D0 = hbm_balance3d('floquet0',hbm,problem,w(i),u(:,i),x(:,i));
    D1 = hbm_balance3d('floquet1',hbm,problem,w(i),u(:,i),x(:,i));
    D2 = hbm_balance3d('floquet2',hbm,problem,w(i),u(:,i),x(:,i));

    I = eye(hbm.harm.NComp*problem.NDof);
    Z = zeros(hbm.harm.NComp*problem.NDof);

    B1 = [D1 D0;
          -I  Z];
    B2 = [D2 Z;
          Z  I];
      
    A(:,:,i) = B1;
    B(:,:,i) = -B2;
end