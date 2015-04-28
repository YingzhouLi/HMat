function A = hscale(alpha, A)
    if A.blockType == 'L'
        A.UMat = sign(alpha)*sqrt(abs(alpha))*A.UMat;
        A.VMat = sqrt(abs(alpha))*A.VMat;
    elseif A.blockType == 'D'
        A.DMat = alpha*A.DMat;
    elseif A.blockType == 'H'
        for i = 1:4
            for j = 1:4
                A.childHMat{i,j} = hscale(alpha,A.childHMat{i,j});
            end
        end
    end
end
