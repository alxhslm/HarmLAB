function X = packdof(Xc,iRetain)
%order is
%   X = [X0;       @ 0
%        real(X1)  @ 1*w0
%        imag(X1)  @ 1*w0
%        real(X2)  @ 2*w0
%        imag(X2)  @ 2*w0
%         ...
%        real(XN)  @ NHarm*w0
%        imag(XN)] @ NHarm*w0

X = mtransposex([real(Xc(2:end,:,:)) imag(Xc(2:end,:,:))]);
X0 = permute(Xc(1,:,:),[2 3 1]);

X = [X0; reshape(X,[],size(X,3))];

if nargin < 2
    iRetain = 1:size(X,1);
end
X = X(iRetain,:);
