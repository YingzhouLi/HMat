function hmatvec(A::HMat2d,v)
    assert(A.width == size(v,1))
    if A.blockType == "LOWRANK"
        return A.UMat*(A.VMat'*v)
    elseif A.blockType == "DENSE"
        return A.DMat*v
    else
        uoffset = 0
        u = zeros(typeof(v[1,1]),A.height,size(v,2))
        for i = 1:4
            voffset = 0
            for j = 1:4
                C = A.childHMat[i,j]
                u[uoffset+(1:C.height),:] += hmatvec(C,v[voffset+(1:C.width),:])
                voffset += C.width
                if j == 4
                    uoffset += C.height
                end
            end
        end
    end
    return u
end

function hmatTvec(A::HMat2d,v)
    assert(A.height == size(v,1))
    if A.blockType == "LOWRANK"
        return A.VMat*(A.UMat'*v)
    elseif A.blockType == "DENSE"
        return A.DMat'*v
    else
        uoffset = 0
        u = zeros(typeof(v[1,1]),A.width,size(v,2))
        for j = 1:4
            voffset = 0
            for i = 1:4
                C = A.childHMat[i,j]
                u[uoffset+(1:C.width),:] += hmatTvec(C,v[voffset+(1:C.height),:])
                voffset += C.height
                if i == 4
                    uoffset += C.width
                end
            end
        end
    end
    return u
end

function hvecmat(v,A::HMat2d)
    assert(A.height == size(v,2))
    if A.blockType == "LOWRANK"
        return v*A.UMat*A.VMat'
    elseif A.blockType == "DENSE"
        return v*A.DMat
    else
        uoffset = 0
        u = zeros(typeof(v[1,1]),size(v,1),A.width)
        for j = 1:4
            voffset = 0
            for i = 1:4
                C = A.childHMat[i,j]
                u[:,uoffset+(1:C.width)] += hvecmat(v[:,voffset+(1:C.height)],C)
                voffset += C.height
                if i == 4
                    uoffset += C.width
                end
            end
        end
    end
    return u
end

function hvecmatT(v,A::HMat2d)
    assert(A.width == size(v,2))
    if A.blockType == "LOWRANK"
        return v*A.VMat*A.UMat'
    elseif A.blockType == "DENSE"
        return v*A.DMat'
    else
        uoffset = 0
        u = zeros(typeof(v[1,1]),size(v,1),A.height)
        for i = 1:4
            voffset = 0
            for j = 1:4
                C = A.childHMat[i,j]
                u[:,uoffset+(1:C.height)] += hvecmatT(v[:,voffset+(1:C.width)],C)
                voffset += C.width
                if j == 4
                    uoffset += C.height
                end
            end
        end
    end
    return u
end
