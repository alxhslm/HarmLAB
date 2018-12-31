function Xc = unpackdof(X,NHarm,NDof,iRetain)
NComp = prod(2*NHarm+1);
if nargin < 4
    iRetain = 1:(NDof*NComp);
end
Xfull = zeros(NDof*NComp,size(X,2));
Xfull(iRetain,:) = X;

if ~isempty(Xfull)
    X = permute(reshape(Xfull,NDof,NComp,[]),[2 1 3]);
else
    X = zeros(NComp,NDof,size(X,2));
end
Xc = [X(1,:,:); X(2:2:end,:,:) + 1i*X(3:2:end,:,:)];