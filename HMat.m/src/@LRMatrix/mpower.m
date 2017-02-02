function C = mpower(A,b)

assert( A.height == A.width );
C = LRMatrix(A.UMat*(A.VMat'*A.UMat)^(b-1),A.VMat,A.EPS,A.MAXRANK);

end