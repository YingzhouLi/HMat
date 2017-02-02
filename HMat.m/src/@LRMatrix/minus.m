function C = minus( A, B )

if isa(A,'LRMatrix') && isa(B,'LRMatrix')
    
    assert( A.height == B.height && A.width == B.width );
    
    C = LRMatrix([A.UMat -B.UMat],[A.VMat B.VMat],A.EPS,A.MAXRANK);
    
else
    C = LR2D(A) - LR2D(B);
end

end