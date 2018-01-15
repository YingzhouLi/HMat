function C = mrdivide( A, B )

if isa(A,'HMatrix') && isa(B,'HMatrix')
    C = mempty(A,B');
    C = hmrdivideh(A,B,C);
elseif isa(A,'LRMatrix') && isa(B,'HMatrix')
    C = LRMatrix(A.UMat,A.VMat'/B,A.EPS,A.MAXRANK);
elseif isa(B,'HMatrix')
    C = dmrdivideh(A,B);
elseif isa(A,'HMatrix')
    C = H2D(A)/B;
end

end

function C = hmrdivideh(A,B,C)

if C.blockType == 'L'
    if B.blockType == 'L'
        error('B cannot be low rank matrix');
    elseif ( B.blockType == 'D' || B.blockType == 'H' ) ...
            && A.blockType == 'L'
        C.LRMat = A.LRMat / B + C.LRMat;
    else
        Rcol = randn(B.height,C.LRMat.MAXRANK+5);
        Rrow = randn(A.height,C.LRMat.MAXRANK+5);
        ABRcol = A*(B\Rcol);
        BTATRrow = B'\(A'*Rrow);
        [Qcol,~] = qr(ABRcol,0);
        [Qrow,~] = qr(BTATRrow,0);
        [Utmp,Stmp,Vtmp] = svdtrunc(...
            pinv(Rrow'*Qcol)*Rrow'*ABRcol*pinv(Qrow'*Rcol),...
            C.LRMat.EPS,C.LRMat.MAXRANK);
        if max(size(Stmp)) > 0
            LRMat = LRMatrix(Qcol*(Utmp*sqrt(Stmp)), ...
                Qrow*(Vtmp*sqrt(Stmp)), C.LRMat.EPS, C.LRMat.MAXRANK );
            C.LRMat = LRMat + C.LRMat;
        end
    end
    C = compress(C);
elseif C.blockType == 'D'
    C.DMat = C.DMat + H2D(A)/H2D(B);
else
    if B.blockType == 'L'
        error('B cannot be low rank matrix');
    elseif ( B.blockType == 'D' || B.blockType == 'H' ) ...
            && A.blockType == 'L'
        C = C + A.LRMat/B;
    elseif B.blockType == 'H' && A.blockType == 'D'
        C = C + A.DMat/B;
    elseif B.blockType == 'D' && A.blockType == 'H'
        C = C + A/B.DMat;
    elseif B.blockType == 'H' && A.blockType == 'H'
        if isupper(B)
            for itj = 1:size(C.childHMat,2)
                for iti = 1:size(C.childHMat,1)
                    Atmp = A.childHMat{iti,itj};
                    for itk = 1:itj-1
                        Atmp = Atmp ...
                            - C.childHMat{iti,itk}*B.childHMat{itk,itj};
                    end
                    C.childHMat{iti,itj} = hmrdivideh( ...
                        Atmp, B.childHMat{itj,itj}, ...
                        C.childHMat{iti,itj});
                end
            end
        elseif islower(B)
            for itj = size(C.childHMat,2):-1:1
                for iti = size(C.childHMat,1):-1:1
                    Atmp = A.childHMat{iti,itj};
                    for itk = size(C.childHMat,2):-1:itj+1
                        Atmp = Atmp ...
                            - C.childHMat{iti,itk}*B.childHMat{itk,itj};
                    end
                    C.childHMat{iti,itj} = hmrdivideh( ...
                        Atmp, B.childHMat{itj,itj}, ...
                        C.childHMat{iti,itj});
                end
            end
        else
            [L,U] = lu(B);
            C = C + (A/U)/L;
        end
    end
end

end

function C = dmrdivideh(A,B)

if B.blockType == 'L'
    error('B cannot be low rank matrix');
elseif B.blockType == 'D'
    C = A/B.DMat;
else
    C = zeros(size(A,1),B.height);
    if isupper(B)
        hoffset = 0;
        for iti = 1:size(B.childHMat,1)
            BC = B.childHMat{iti,iti};
            Atmp = A(:,hoffset+(1:BC.width));
            offset = 0;
            for itj = 1:iti-1
                Cidx = offset + (1:B.childHMat{itj,iti}.height);
                Atmp = Atmp ...
                    - C(:,Cidx)*B.childHMat{itj,iti};
                offset = Cidx(end);
            end
            C(:,hoffset+(1:BC.height)) = Atmp / BC;
            hoffset = hoffset + BC.width;
        end
    elseif islower(B)
        hoffset = B.height+1;
        for iti = size(B.childHMat,1):-1:1
            BC = B.childHMat{iti,iti};
            Atmp = A(:,hoffset-(BC.width:-1:1));
            offset = B.height+1;
            for itj = size(B.childHMat,1):-1:iti+1
                Cidx = offset - (B.childHMat{itj,iti}.height:-1:1);
                Atmp = Atmp ...
                    - C(:,Cidx)*B.childHMat{itj,iti};
                offset = Cidx(1);
            end
            C(:,hoffset-(BC.height:-1:1),:) = Atmp/ BC;
            hoffset = hoffset - BC.width;
        end
    else
        [L,U] = lu(B);
        C = C + (A/U)/L;
    end
end

end
