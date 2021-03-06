function +(A::LRMat,B::LRMat)
    assert(A.height == B.height);
    assert(A.width == B.width);

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
        if B.blockType == LOWRANK
            C.UMat = [A.UMat B.UMat];
            C.VMat = [A.VMat B.VMat];
        else
            (Utmp,Vtmp) = hh2l(B);
            C.UMat = [A.UMat Utmp];
            C.VMat = [A.VMat Vtmp];
        end
        hcompress(C);
    elseif A.blockType == DENSE
        C.blockType = DENSE;
        C.DMat = A.DMat + hh2d(B);
    elseif A.blockType == HMAT
        C.blockType = HMAT;
        C.childHMat = Array(HMat{T},4,4);
        hoffset = 0;
        for i = 1:4
            woffset = 0;
            for j = 1:4
                if B.blockType == LOWRANK
                    BH = HMat();
                    BH.height = A.childHMat[i,j].height;
                    BH.width = A.childHMat[i,j].width;
                    BH.blockType = LOWRANK;
                    BH.UMat = B.UMat[hoffset+(1:BH.height),:];
                    BH.VMat = B.VMat[woffset+(1:BH.width),:];
                    woffset += BH.width;
                    if j == 4
                        hoffset += BH.height;
                    end
                    C.childHMat[i,j] = hadd(A.childHMat[i,j],BH);
                elseif B.blockType == DENSE
                    BH = HMat();
                    BH.height = A.childHMat[i,j].height;
                    BH.width = A.childHMat[i,j].width;
                    BH.blockType = DENSE;
                    BH.EPS = B.EPS;
                    BH.MAXRANK = B.MAXRANK;
                    BH.MINN = B.MINN;
                    BH.DMat = B.DMat[hoffset+(1:BH.height),woffset+(1:BH.width)];
                    woffset += BH.width;
                    if j == 4
                        hoffset += BH.height;
                    end
                    C.childHMat[i,j] = hadd(A.childHMat[i,j],BH);
                else
                    C.childHMat[i,j] = hadd(A.childHMat[i,j],B.childHMat[i,j]);
                end
            end
        end
    end
    return C;
end

function hadd{T<:Number}(A::HMat{T},UMat::Matrix{T},VMat::Matrix{T})
    assert(A.height == size(UMat,1));
    assert(A.width == size(VMat,1));

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
        C.UMat = [A.UMat UMat];
        C.VMat = [A.VMat VMat];
        hcompress(C);
    elseif A.blockType == DENSE
        C.blockType = DENSE;
        C.DMat = A.DMat + UMat*VMat';
    elseif A.blockType == HMAT
        C.blockType = HMAT;
        C.childHMat = Array(HMat{T},4,4);
        hoffset = 0;
        for i = 1:4
            woffset = 0;
            for j = 1:4
                C.childHMat[i,j] = hadd(A.childHMat[i,j],UMat[hoffset+(1:A.childHMat[i,j].height),:],VMat[woffset+(1:A.childHMat[i,j].width),:]);
                woffset += A.childHMat[i,j].width;
                if j == 4
                    hoffset += A.childHMat[i,j].height;
                end
            end
        end
    end
    return C;
end

function hadd!(A::HMat,UMat::AbstractMatrix,VMat::AbstractMatrix)
    assert(A.height == size(UMat,1));
    assert(A.width == size(VMat,1));

    if A.blockType == LOWRANK
        A.UMat = [A.UMat UMat];
        A.VMat = [A.VMat VMat];
        hcompress(A);
    elseif A.blockType == DENSE
        A.DMat += UMat*VMat';
    elseif A.blockType == HMAT
        hoffset = 0;
        for i = 1:4
            woffset = 0;
            for j = 1:4
                urange = hoffset .+ (1:A.childHMat[i,j].height)
                vrange = woffset .+ (1:A.childHMat[i,j].width)
                hadd!(A.childHMat[i,j],
                      sub(UMat,urange,:),
                      sub(VMat,vrange,:))
                woffset += A.childHMat[i,j].width;
                if j == 4
                    hoffset += A.childHMat[i,j].height;
                end
            end
        end
    end
end

function hadddiag{T<:Number}(A::HMat{T},alpha)
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
        C.UMat = A.UMat;
        C.VMat = A.VMat;
    elseif A.blockType == DENSE
        C.blockType = DENSE;
        C.DMat = A.DMat;
        if C.trg == C.src
            C.DMat += alpha*eye(C.height);
        end
    elseif A.blockType == HMAT
        C.blockType = HMAT;
        C.childHMat = Array(HMat{T},4,4);
        for i = 1:4, j = 1:4
            C.childHMat[i,j] = hadddiag(A.childHMat[i,j],alpha);
        end
    end
    return C;
end

function hadddiag!(A::HMat,alpha)
    if A.blockType == DENSE
        if A.trg == A.src
            A.DMat += alpha*eye(A.height);
        end
    elseif A.blockType == HMAT
        for i = 1:4
            hadddiag!(A.childHMat[i,i],alpha);
        end
    end
end
