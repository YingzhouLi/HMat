function C = mpower(A,b)

C = LRMatrix(A.UMat*(A.VMat'*A.UMat)^(b-1),A.VMat,A.EPS,A.MAXRANK);

end