function y = normfro(A)

if A.blockType == 'L'
    LRMat = A.LRMat;
    y = sqrt(trace(LRMat.UMat'*LRMat.UMat*(LRMat.VMat'*LRMat.VMat)));
elseif A.blockType == 'D'
    y = norm(A.DMat,'fro');
elseif A.blockType == 'H'
    y = 0;
    for iti = 1:size(A.childHMat,1)
        for itj = 1:size(A.childHMat,2)
            y = y + normfro(A.childHMat{iti,itj})^2;
        end
    end
    y = sqrt(y);
end

end