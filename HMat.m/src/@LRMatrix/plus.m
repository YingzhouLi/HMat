function C = plus( A, B )

if isa(A,'LRMatrix') && isa(B,'LRMatrix')
    C = LRMatrix([A.UMat B.UMat],[A.VMat B.VMat],A.EPS,A.MAXRANK);
elseif isa(A,'LRMatrix')
    C = LR2D(A) + B;
elseif isa(B,'LRMatrix')
    C = A + LR2D(B);
end

end