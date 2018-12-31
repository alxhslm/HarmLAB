function C = blkmat(A)
B = num2cell(A,[1,2]);
C = blkdiag(B{:});