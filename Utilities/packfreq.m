function X = packfreq(Xfull,iSub)
X = zeros(length(iSub),size(Xfull,3));

for i = 1:length(iSub)
    X(i,:) = Xfull(iSub(i,1),iSub(i,2),:);
end

%now double everything as we want single-sided coefficients
X(2:end,:) = 2*X(2:end,:); %except 0Hz
