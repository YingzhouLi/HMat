function C = empty(A)

C = A;
if A.blockType == 'L'
    C.LRMat = LRMatrix( zeros(A.height,0), zeros(A.width,0), ...
        A.LRMat.EPS, A.LRMat.MAXRANK );
elseif A.blockType == 'D'
    C.DMat = zeros(A.height,A.width);
elseif A.blockType == 'H'
    C.childHMat = cell(size(A.childHMat));
    for i = 1:length(A.childHMat)
        C.childHMat{i} = empty(A.childHMat{i});
    end
end

end