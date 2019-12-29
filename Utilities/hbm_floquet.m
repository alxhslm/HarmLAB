function lambda = hbm_floquet(hbm,problem,w,U,X)
if ~(isvector(X) && size(X,1) == hbm.harm.NRetainNL)
    x = packdof(X(:,problem.iNL),hbm.harm.iRetainNL);
    u = packdof(U);
else
    x = X;
    u = U;
end

if any(isnan(x) | isinf(x))
    lambda = NaN(hbm.harm.NRetainNL,1);
    return
end

[A,B] = floquetMatrices(hbm,problem,w,u,x);
lambda = eig(A,B,'vector');
[~,iSort] = sort(abs(imag(lambda)));
lambda = lambda(iSort);

% %remove real eigenvalues
% iKeep = abs(imag(lambda)) > 1E-6;
% lambda = lambda(iKeep);

% lambda = lambda((1:2*problem.NDof),:);

function [A,B] = floquetMatrices(hbm,problem,w,u,x)
hbm.bIncludeNL = 1;
NPts = size(u,2);
% h = waitbar(0,'Computing matrices');

A = zeros(2*size(x,1),2*size(x,1),NPts);
B = zeros(2*size(x,1),2*size(x,1),NPts);
for i = 1:NPts
    D0 = hbm_balance3d('floquet0',hbm,problem,w(i),u(:,i),x(:,i));
    D1 = hbm_balance3d('floquet1',hbm,problem,w(i),u(:,i),x(:,i));
    D2 = hbm_balance3d('floquet2',hbm,problem,w(i),u(:,i),x(:,i));

    % lambda = polyeig(D0,D1,D2);

    I = eye(hbm.harm.NRetainNL);
    Z = zeros(hbm.harm.NRetainNL);

    B1 = [D1 D0;
          -I  Z];
    B2 = [D2 Z;
          Z  I];
      
    A(:,:,i) = B1;
    B(:,:,i) = -B2;
    
%     waitbar(i/NPts,h)
end
% close(h)