function hidentity!(A::HMat)
    if A.blockType == LOWRANK
        A.UMat = zeros(A.height,0);
        A.VMat = zeros(A.width,0);
    elseif A.blockType == DENSE
        if A.trg == A.src
            A.DMat = eye(A.height,A.width)
        else
            A.DMat = zeros(A.height,A.width);
        end
    elseif A.blockType == HMAT
        for i = 1:4, j = 1:4
            hidentity!(A.childHMat[i,j]);
        end
    end
end

function hidentity{T<:Number}(A::HMat{T})
    C = HMat{T}();
    C.height = A.height;
    C.width = A.width;
    C.trg = A.trg;
    C.src = A.src;
    C.level = A.level;

    C.EPS = A.EPS;
    C.MAXRANK = A.MAXRANK;
    C.MINN = A.MINN;

    if A.blockType == LOWRANK
        C.blockType = LOWRANK;
        C.UMat = zeros(C.height,0);
        C.VMat = zeros(C.width,0);
    elseif A.blockType == DENSE
        C.blockType = DENSE;
        if C.trg == C.src
            C.DMat = eye(C.height,C.width);
        else
            C.DMat = zeros(C.height,C.width);
        end
    elseif A.blockType == HMAT
        C.blockType = HMAT;
        C.childHMat = Array(HMat{T},4,4);
        for i = 1:4
            for j = 1:4
                C.childHMat[i,j] = hidentity(A.childHMat[i,j]);
            end
        end
    end
    return C;
end
