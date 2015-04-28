function C = hmuladd(A,B,C)
    assert(A.width == B.height);
    assert(A.height == C.height);
    assert(B.width == C.width);

    if C.blockType == 'L'
        if A.blockType == 'L' && B.blockType == 'L'
            Mtmp = A.VMat'*B.UMat;
            [Utmp,Stmp,Vtmp] = svdtrunc(Mtmp,C.EPS,C.MAXRANK);
            if max(size(Stmp)) > 0
                C.UMat = [C.UMat A.UMat*(Utmp*sqrt(Stmp))];
                C.VMat = [C.VMat B.VMat*(Vtmp*sqrt(Stmp))];
            end
        elseif A.blockType == 'L' && ( B.blockType == 'D' || B.blockType == 'H' )
            C.UMat = [C.UMat A.UMat];
            C.VMat = [C.VMat hmatTvec(B,A.VMat)];
        elseif ( A.blockType == 'D' || A.blockType == 'H' ) && B.blockType == 'L'
            C.UMat = [C.UMat hmatvec(A,B.UMat)];
            C.VMat = [C.VMat B.VMat];
        else
            Rcol = randn(B.width,C.MAXRANK+5);
            Rrow = randn(A.height,C.MAXRANK+5);
            ABRcol = hmatvec(A,hmatvec(B,Rcol));
            BTATRrow = hmatTvec(B,hmatTvec(A,Rrow));
            [Qcol,~] = qr(ABRcol,0);
            [Qrow,~] = qr(BTATRrow,0);
            [Utmp,Stmp,Vtmp] = svdtrunc(...
                        pinv(Rrow'*Qcol)*Rrow'*ABRcol*pinv(Qrow'*Rcol),...
                        C.EPS,C.MAXRANK);
            if max(size(Stmp)) > 0
                C.UMat = [C.UMat Qcol*(Utmp*sqrt(Stmp))];
                C.VMat = [C.VMat Qrow*(Vtmp*sqrt(Stmp))];
            end
        end
        C = hcompress(C);
    elseif C.blockType == 'D'
        C.DMat = C.DMat + hh2d(A)*hh2d(B);
    else
        if A.blockType == 'L' && B.blockType == 'L'
            Mtmp = A.VMat'*B.UMat;
            [Utmp,Stmp,Vtmp] = svdtrunc(Mtmp,C.EPS,C.MAXRANK);
            if max(size(Stmp)) > 0
                C = haddl(C,A.UMat*(Utmp*sqrt(Stmp)),...
                            B.VMat*(Vtmp*sqrt(Stmp)));
            end
        elseif A.blockType == 'L' && ( B.blockType == 'D' || B.blockType == 'H' )
            C = haddl(C,A.UMat,hmatTvec(B,A.VMat));
        elseif ( A.blockType == 'D' || A.blockType == 'H' ) && B.blockType == 'L'
            C = haddl(C,hmatvec(A,B.UMat),B.VMat);
        elseif A.blockType == 'H' && B.blockType == 'D'
            error('Not been implemented');
        elseif A.blockType == 'D' && B.blockType == 'H'
            error('Not been implemented');
        elseif A.blockType == 'H' && B.blockType == 'H'
            for i = 1:4
            for j = 1:4
            for k = 1:4
                C.childHMat{i,j} = hmuladd(A.childHMat{i,k},B.childHMat{k,j},...
                            C.childHMat{i,j});
            end
            end
            end
        end
    end
end
