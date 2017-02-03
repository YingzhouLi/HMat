function C = mtimes(A,B)

if isa(A,'LRMatrix') && isa(B,'LRMatrix')
    
    if size(A.UMat,2) <= size(B.UMat,2)
        C = LRMatrix(A.UMat,B.VMat*(B.UMat'*A.VMat),A.EPS,A.MAXRANK);
    else
        C = LRMatrix(A.UMat*(A.VMat'*B.UMat),B.VMat,A.EPS,A.MAXRANK);
    end
    
elseif isa(A,'LRMatrix')
    C = LRMatrix(A.UMat,B'*A.VMat,A.EPS,A.MAXRANK);
elseif isa(B,'LRMatrix')
    C = LRMatrix(A*B.UMat,B.VMat,B.EPS,B.MAXRANK);
end

end