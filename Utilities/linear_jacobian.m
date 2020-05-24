function [Bk,Bc,Bm,Bx] = linear_jacobian(harm,Ku,Cu,Mu,Gu)
if nargin < 5
    Gu = 0*Cu;
end

%terms dependent on each frequency
Ak{1} = Ku;
for j = 2:harm.NFreq
    Ak{j} = blkdiag(Ku,Ku);
end

for k = 1:2
    Ac{k}{1} = 0*Cu;
    Am{k}{1} = 0*Mu;
    for j = 2:harm.NFreq
        Ac{k}{j} =  harm.rFreqBase(k)*harm.kHarm(j,k)*antiblkdiag(-Cu,Cu);
        Am{k}{j} = -(harm.rFreqBase(k)*harm.kHarm(j,k))^2*blkdiag(Mu,Mu);
    end
end

%cross-terms
Ax{1} = 0*Mu;
for j = 2:harm.NFreq
    Ax{j} = -2*prod(harm.kHarm(j,:).*harm.rFreqBase)*blkdiag(Mu,Mu);
end

%add on gyro terms
for j = 2:harm.NFreq
    Am{1}{j} = Am{1}{j} + harm.rFreqBase(1)*harm.kHarm(j,1)*antiblkdiag(-Gu,Gu);
    Ax{j}    = Ax{j}    + harm.rFreqBase(2)*harm.kHarm(j,2)*antiblkdiag(-Gu,Gu);
end

Bk = blkdiag(Ak{:});
for k = 1:2
    Bc{k} = blkdiag(Ac{k}{:});
    Bm{k} = blkdiag(Am{k}{:});
end
Bx = blkdiag(Ax{:});

Bk = sparse(Bk);
for k = 1:2
    Bm{k} = sparse(Bm{k});
    Bc{k} = sparse(Bc{k});
end
Bx = sparse(Bx);
