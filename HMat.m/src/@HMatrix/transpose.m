function C = transpose( A )

C = A;
C.height = A.width;
C.width = A.height;
if A.blockType == 'L'
    C.LRMat = A.LRMat.';
elseif A.blockType == 'D'
    C.DMat = A.DMat.';
elseif A.blockType == 'H'
    C.childHMat = cell(size(A.childHMat,2),size(A.childHMat,1));
    for iti = 1:size(A.childHMat,1)
        for itj = 1:size(A.childHMat,2)
            C.childHMat{itj,iti} = A.childHMat{iti,itj}.';
        end
    end
end

end