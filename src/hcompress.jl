function hcompress(A::HMat2d)
    if A.blockType == "LOWRANK"
        (QU,RU) = qr(A.UMat)
        (QV,RV) = qr(A.VMat)
        (Utmp,Stmp,Vtmp) = svdtrunc(RU*RV',A.EPS,A.MAXRANK)
        A.rank = length(Stmp)
        A.UMat = QU*(Utmp.*sqrt(Stmp)')
        A.VMat = QV*(Vtmp.*sqrt(Stmp)')
    elseif A.blockType == "HMAT"
        for C in A.childHMat
            hcompress(C)
        end
    end
end

function hh2l(A::HMat2d)
    Rcol = randn(A.width,A.MAXRANK+5)
    Rrow = randn(A.height,A.MAXRANK+5)
    ARcol = hmatvec(A,Rcol)
    ATRrow = hmatTvec(A,Rrow)
    (Qcol,) = qr(ARcol)
    (Qrow,) = qr(ATRrow)
    (Utmp,Stmp,Vtmp) = svdtrunc(pinv(Rrow'*Qcol)*Rrow'*ARcol*pinv(Qrow'*Rcol),A.MAXRANK,A.EPS)
    UMat = Qcol*(Utmp.*sqrt(Stmp)')
    VMat = Qrow*(Vtmp.*sqrt(Stmp)')
    return UMat,VMat
end
