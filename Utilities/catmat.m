function C = catmat(A,dim)
B = num2cell(A,[1,2]);
C = cat(dim,B{:});