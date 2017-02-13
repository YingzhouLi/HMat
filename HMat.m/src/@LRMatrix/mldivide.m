function C = mldivide( A, B )

if isa(A,'LRMatrix')
    error('LR matrix cannot be the denominator');
end

C = LRMatrix(A\B.UMat,B.VMat,B.EPS,B.MAXRANK);

end