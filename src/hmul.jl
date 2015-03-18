function hmul(A::HMat2d,B::HMat2d,C::HMat2d)
    assert(A.width == B.height)
    assert(A.height == C.height)
    assert(B.width == C.width)
    assert(C.trg == A.trg)
    assert(C.src == B.src)

    if C.blockType == LOWRANK
        if A.blockType == LOWRANK && B.blockType == LOWRANK
            Mtmp = A.VMat'*B.UMat
            (Utmp,Stmp,Vtmp) = svdtrunc(Mtmp,C.EPS,C.MAXRANK)
            if length(Stmp) > 0
                C.UMat = [C.UMat A.UMat*(Utmp.*sqrt(Stmp)')]
                C.VMat = [C.VMat B.VMat*(Vtmp.*sqrt(Stmp)')]
            end
        elseif A.blockType == LOWRANK && ( B.blockType == DENSE || B.blockType == HMAT )
            C.UMat = [C.UMat A.UMat]
            C.VMat = [C.VMat hmatTvec(B,A.VMat)]
        elseif ( A.blockType == DENSE || A.blockType == HMAT ) && B.blockType == LOWRANK
            C.UMat = [C.UMat hmatvec(A,B.UMat)]
            C.VMat = [C.VMat B.VMat]
        else
            Rcol = randn(B.width,C.MAXRANK+5)
            Rrow = randn(A.height,C.MAXRANK+5)
            ABRcol = hmatvec(A,hmatvec(B,Rcol))
            BTATRrow = hmatTvec(B,hmatTvec(A,Rrow))
            (Qcol,) = qr(ABRcol)
            (Qrow,) = qr(BTATRrow)
            (Utmp,Stmp,Vtmp) = svdtrunc(pinv(Rrow'*Qcol)*Rrow'*ABRcol*pinv(Qrow'*Rcol),C.EPS,C.MAXRANK)
            if length(Stmp) > 0
                C.UMat = [C.UMat Qcol*(Utmp.*sqrt(Stmp)')]
                C.VMat = [C.VMat Qrow*(Vtmp.*sqrt(Stmp)')]
            end
        end
        hcompress(C)
    elseif C.blockType == DENSE
        C.DMat += hh2d(A)*hh2d(B)
    else
        if A.blockType == LOWRANK && B.blockType == LOWRANK
            Mtmp = A.VMat'*B.UMat
            (Utmp,Stmp,Vtmp) = svdtrunc(Mtmp,C.EPS,C.MAXRANK)
            if length(Stmp) > 0
                hadd!(C,A.UMat*(Utmp.*sqrt(Stmp)'),B.VMat*(Vtmp.*sqrt(Stmp)'))
            end
        elseif A.blockType == LOWRANK && ( B.blockType == DENSE || B.blockType == HMAT )
            hadd!(C,A.UMat,hmatTvec(B,A.VMat))
        elseif ( A.blockType == DENSE || A.blockType == HMAT ) && B.blockType == LOWRANK
            hadd!(C,hmatvec(A,B.UMat),B.VMat)
        elseif A.blockType == HMAT && B.blockType == DENSE
            error("Not been implemented")
        elseif A.blockType == DENSE && B.blockType == HMAT
            error("Not been implemented")
        elseif A.blockType == HMAT && B.blockType == HMAT
            for i = 1:4, j = 1:4, k = 1:4
                hmul(A.childHMat[i,k],B.childHMat[k,j],C.childHMat[i,j])
            end
        end
    end
end

function hmul(A::HMat2d,B::HMat2d)
    assert(A.width == B.height)
    C = hempty(A)
    hmul(A,B,C)
    return C
end
