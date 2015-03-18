function hadjoint(A::HMat2d)
    C = HMat2d();
    C.height = A.width;
    C.width = A.height;
    C.trg = A.src;
    C.src = A.trg;
    C.level = A.level;

    C.EPS = A.EPS;
    C.MAXRANK = A.MAXRANK;
    C.MINN = A.MINN;

    if A.blockType == LOWRANK
        C.blockType = LOWRANK;
        C.UMat = A.VMat;
        C.VMat = A.UMat;
    elseif A.blockType == DENSE
        C.blockType = DENSE;
        C.DMat = A.DMat';
    elseif A.blockType == HMAT
        C.blockType = HMAT;
        C.childHMat = Array(HMat2d,4,4);
        for i = 1:4
            for j = 1:4
                C.childHMat[j,i] = hadjoint(A.childHMat[i,j]);
            end
        end
    end
    return C;
end
