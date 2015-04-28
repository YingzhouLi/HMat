function y = hnormfro(A)
    if A.blockType == 'L'
        y = sqrt(trace(A.UMat'*A.UMat*(A.VMat'*A.VMat)));
    elseif A.blockType == 'D'
        y = norm(A.DMat,'fro');
    elseif A.blockType == 'H'
        y = 0;
        for i = 1:4
            for j = 1:4
                y = y + hnormfro(A.childHMat{i,j})^2;
            end
        end
        y = sqrt(y);
    end
end
