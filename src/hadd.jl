function hadd(A::HMat2d,B::HMat2d)
    assert(A.height == B.height)
    assert(A.width == B.width)

    C = HMat2d()
    C.height = A.height
    C.width = A.width
    C.trg = A.trg
    C.src = A.src
    C.level = A.level

    C.EPS = A.EPS
    C.MAXRANK = A.MAXRANK
    C.MINN = A.MINN

    if A.blockType == "LOWRANK"
        C.blockType = "LOWRANK"
        if B.blockType == "LOWRANK"
            C.UMat = [A.UMat B.UMat]
            C.VMat = [A.VMat B.VMat]
        else
            (Utmp,Vtmp) = hh2l(B)
            C.UMat = [A.UMat Utmp]
            C.VMat = [A.VMat Vtmp]
        end
        hcompress(C)
    elseif A.blockType == "DENSE"
        C.blockType = "DENSE"
        C.DMat = A.DMat + hh2d(B)
    elseif A.blockType == "HMAT"
        C.blockType = "HMAT"
        C.childHMat = Array(HMat2d,4,4)
        hoffset = 0
        for i = 1:4
            woffset = 0
            for j = 1:4
                if B.blockType == "LOWRANK"
                    BH = HMat2d()
                    BH.height = A.childHMat[i,j].height
                    BH.width = A.childHMat[i,j].width
                    BH.blockType = "LOWRANK"
                    BH.UMat = B.UMat[hoffset+(1:BH.height),:]
                    BH.VMat = B.VMat[woffset+(1:BH.width),:]
                    woffset += BH.width
                    if j == 4
                        hoffset += BH.height
                    end
                    C.childHMat[i,j] = hadd(A.childHMat[i,j],BH)
                elseif B.blockType == "DENSE"
                    BH = HMat2d()
                    BH.height = A.childHMat[i,j].height
                    BH.width = A.childHMat[i,j].width
                    BH.blockType = "DENSE"
                    BH.EPS = B.EPS
                    BH.MAXRANK = B.MAXRANK
                    BH.MINN = B.MINN
                    BH.DMat = B.DMat[hoffset+(1:BH.height),woffset+(1:BH.width)]
                    woffset += BH.width
                    if j == 4
                        hoffset += BH.height
                    end
                    C.childHMat[i,j] = hadd(A.childHMat[i,j],BH)
                else
                    C.childHMat[i,j] = hadd(A.childHMat[i,j],B.childHMat[i,j])
                end
            end
        end
    end
    return C
end

function hadd(A::HMat2d,UMat::Array,VMat::Array)
    assert(A.height == size(UMat,1))
    assert(A.width == size(VMat,1))

    C = HMat2d()
    C.height = A.height
    C.width = A.width
    C.trg = A.trg
    C.src = A.src
    C.level = A.level

    C.EPS = A.EPS
    C.MAXRANK = A.MAXRANK
    C.MINN = A.MINN

    if A.blockType == "LOWRANK"
        C.blockType = "LOWRANK"
        C.UMat = [A.UMat UMat]
        C.VMat = [A.VMat VMat]
        hcompress(C)
    elseif A.blockType == "DENSE"
        C.blockType = "DENSE"
        C.DMat = A.DMat + UMat*VMat'
    elseif A.blockType == "HMAT"
        C.blockType = "HMAT"
        C.childHMat = Array(HMat2d,4,4)
        hoffset = 0
        for i = 1:4
            woffset = 0
            for j = 1:4
                C.childHMat[i,j] = hadd(A.childHMat[i,j],UMat[hoffset+(1:A.childHMat[i,j].height),:],VMat[woffset+(1:A.childHMat[i,j].width),:])
                woffset += A.childHMat[i,j].width
                if j == 4
                    hoffset += A.childHMat[i,j].height
                end
            end
        end
    end
    return C
end

function hadd!(A::HMat2d,UMat::Array,VMat::Array)
    assert(A.height == size(UMat,1))
    assert(A.width == size(VMat,1))

    if A.blockType == "LOWRANK"
        A.UMat = [A.UMat UMat]
        A.VMat = [A.VMat VMat]
        hcompress(A)
    elseif A.blockType == "DENSE"
        A.DMat += UMat*VMat'
    elseif A.blockType == "HMAT"
        hoffset = 0
        for i = 1:4
            woffset = 0
            for j = 1:4
                hadd!(A.childHMat[i,j],UMat[hoffset+(1:A.childHMat[i,j].height),:],VMat[woffset+(1:A.childHMat[i,j].width),:])
                woffset += A.childHMat[i,j].width
                if j == 4
                    hoffset += A.childHMat[i,j].height
                end
            end
        end
    end
end
