function D = hh2d(A)
    if A.blockType == 'D'
        D = A.DMat;
    elseif A.blockType == 'L'
        D = A.UMat*A.VMat';
    elseif A.blockType == 'H'
        D = zeros(A.height,A.width);
        hoffset = 0;
        for i = 1:4
            woffset = 0;
            for j = 1:4
                C = A.childHMat{i,j};
                D(hoffset+(1:C.height),woffset+(1:C.width)) = hh2d(C);
                woffset = woffset + C.width;
                if j == 4
                    hoffset = hoffset + C.height;
                end
            end
        end
    end
end
