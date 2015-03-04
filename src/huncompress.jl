function hh2d(A::HMat2d)
    if A.blockType == "LOWRANK"
        return A.UMat*A.VMat'
    elseif A.blockType == "DENSE"
        return A.DMat
    elseif A.blockType == "HMAT"
        D = zeros(A.height,A.width)
        hoffset = 0
        for i = 1:4
            woffset = 0
            for j = 1:4
                C = A.childHMat[i,j]
                D[hoffset+(1:C.height),woffset+(1:C.width)] = hh2d(C)
                woffset += C.width
                if j == 4
                    hoffset += C.height
                end
            end
        end
        return D
    end
end
