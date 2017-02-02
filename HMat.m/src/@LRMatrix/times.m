function C = times(A,B)

if isa(A,'LRMatrix') && isa(B,'LRMatrix')
    
    assert( A.height == B.height && A.width == B.width );
    
    rA = size(A.UMat,2);
    rB = size(B.UMat,2);
    UMat = reshape( repmat(A.UMat,rB,1), A.height,[] ) ...
        .* repmat(B.UMat,1,rA);
    VMat = reshape( repmat(A.VMat,rB,1), A.width,[] ) ...
        .* repmat(B.VMat,1,rA);
    C = LRMatrix(UMat,VMat,A.EPS,A.MAXRANK);
    
else
    C = LR2D(A).*LR2D(B);
end

end