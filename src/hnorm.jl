function hnorm(A::HMat2d,p="fro")
    if p == "fro"
        return hnormfro(A)
    elseif p == Inf
        error("Not been implemented")
    elseif p == 1
        error("Not been implemented")
    elseif p == 2
        error("Not been implemented")
    else
        error("hnorm only support p = 1,2,Inf")
    end
end

function hnormfro(A::HMat2d)
    if A.blockType == LOWRANK
        return sqrt(trace(A.UMat'*A.UMat*(A.VMat'*A.VMat)))
    elseif A.blockType == DENSE
        return vecnorm(A.DMat)
    elseif A.blockType == HMAT
        fro = 0
        for i = 1:4
            for j = 1:4
                fro += hnormfro(A.childHMat[i,j])^2
            end
        end
        return sqrt(fro)
    end
end
