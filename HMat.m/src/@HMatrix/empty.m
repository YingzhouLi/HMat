function C = empty(A)

C = A;
if A.blockType == 'L'
    C.LRMat = LRMatrix( zeros(A.height,0), zeros(A.width,0), ...
        A.LRMat.EPS, A.LRMat.MAXRANK );
elseif A.blockType == 'D'
    C.DMat = zeros(A.height,A.width);
elseif A.blockType == 'H'
    for it = 1:numel(A.childHMat)
        C.childHMat{it} = empty(A.childHMat{it});
    end
else
    error('Wrong');
end

end