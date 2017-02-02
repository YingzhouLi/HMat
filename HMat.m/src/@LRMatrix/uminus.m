function C = uminus( A )

C = LRMatrix(-A.UMat,A.VMat,A.EPS,A.MAXRANK);

end