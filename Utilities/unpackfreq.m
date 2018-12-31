function Xfull = unpackfreq(X,NHarm,iSub)
Xfull = zeros(2*NHarm(1)+1,2*NHarm(2)+1,size(X,2));
for i = 1:length(iSub)
    Xfull(iSub(i,1),iSub(i,2),:) = X(i,:);
end
Xfull = 0.5*(Xfull + conj(flip(flip(Xfull,1),2)));