function C = uminus( A )

C = A;
if A.blockType == 'L'
    C.LRMat = - A.LRMat;
elseif A.blockType == 'D'
    C.DMat = A.DMat;
elseif A.blockType == 'H'
    for iti = 1:size(A.childHMat,1)
        for itj = 1:size(A.childHMat,2)
            C.childHMat{iti,itj} = - A.childHMat{iti,itj};
        end
    end
end

end