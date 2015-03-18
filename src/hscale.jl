function hscale(alpha, A::HMat2d)
    C = HMat2d()
    C.height = A.height
    C.width = A.width
    C.trg = A.trg
    C.src = A.src
    C.level = A.level
    C.blockType = A.blockType

    C.EPS = A.EPS
    C.MAXRANK = A.MAXRANK
    C.MINN = A.MINN

    if A.blockType == LOWRANK
        C.UMat = sign(alpha)*sqrt(abs(alpha))*A.UMat
        C.VMat = sqrt(abs(alpha))*A.VMat
    elseif A.blockType == DENSE
        C.DMat = alpha*A.DMat
    elseif A.blockType == HMAT
        C.blockType = HMAT
        C.childHMat = Array(HMat2d,4,4)
        for i = 1:4, j = 1:4
            C.childHMat[i,j] = hscale(alpha,A.childHMat[i,j])
        end
    end
    return C
end

function hscale!(alpha, A::HMat2d)
    if A.blockType == LOWRANK
        A.UMat *= sign(alpha)*sqrt(abs(alpha))
        A.VMat *= sqrt(abs(alpha))
    elseif A.blockType == DENSE
        A.DMat *= alpha
    elseif A.blockType == HMAT
        for C in A.childHMat
            hscale!(alpha,C)
        end
    end
end

