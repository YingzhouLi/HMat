function D = H2D(A)

if A.blockType == 'D'
    D = A.DMat;
elseif A.blockType == 'L'
    D = LR2D(A.LRMat);
elseif A.blockType == 'H'
    D = zeros(A.height,A.width);
    hoffset = 0;
    for iti = 1:size(A.childHMat,1)
        woffset = 0;
        for itj = 1:size(A.childHMat,2)
            C = A.childHMat{iti,itj};
            D(hoffset+(1:C.height),woffset+(1:C.width)) = H2D(C);
            woffset = woffset + C.width;
        end
        hoffset = hoffset + C.height;
    end
else
    D = A;
end

end
