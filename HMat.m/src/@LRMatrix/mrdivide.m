function C = mrdivide( A, B )

if isa(B,'LRMatrix')
    error('LR matrix cannot be the denominator');
end

C = LRMatrix(A.UMat,B'\A.VMat,A.EPS,A.MAXRANK);

end