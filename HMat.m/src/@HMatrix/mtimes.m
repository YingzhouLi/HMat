function C = mtimes( A, B )

if isa(A,'HMatrix') && isa(B,'HMatrix')
    C = empty(A);
    C = hmtimesh(A,B,C);
elseif isa(A,'HMatrix') && isa(B,'LRMatrix')
    C = LRMatrix(A*B.UMat,B.VMat,B.EPS,B.MAXRANK);
elseif isa(A,'LRMatrix') && isa(B,'HMatrix')
    C = LRMatrix(A.UMat,A.VMat'*B,A.EPS,A.MAXRANK);
elseif isa(A,'HMatrix')
    C = hmtimesd(A,B);
elseif isa(B,'HMatrix')
    C = dmtimesh(A,B);
end

end

function C = hmtimesh(A,B,C)

if A.blockType == 'E' || B.blockType == 'E'
    return;
end

if C.blockType == 'L'
    if A.blockType == 'L' && B.blockType == 'L'
        C.LRMat = A.LRMat*B.LRMat + C.LRMat;
    elseif A.blockType == 'L' ...
            && ( B.blockType == 'D' || B.blockType == 'H' )
        C.LRMat = A.LRMat*B + C.LRMat;
    elseif ( A.blockType == 'D' || A.blockType == 'H' ) ...
            && B.blockType == 'L'
        C.LRMat = A*B.LRMat + C.LRMat;
    else
        Rcol = randn(B.width,C.LRMat.MAXRANK+5);
        Rrow = randn(A.height,C.LRMat.MAXRANK+5);
        ABRcol = A*(B*Rcol);
        BTATRrow = B'*(A'*Rrow);
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
    C.DMat = C.DMat + H2D(A)*H2D(B);
else
    if A.blockType == 'L' && B.blockType == 'L'
        C = C + A.LRMat*B.LRMat;
    elseif A.blockType == 'L' && ...
            ( B.blockType == 'D' || B.blockType == 'H' )
        C = C + A.LRMat*B;
    elseif ( A.blockType == 'D' || A.blockType == 'H' ) ...
            && B.blockType == 'L'
        C = C + A*B.LRMat;
    elseif A.blockType == 'H' && B.blockType == 'D'
        error('Not been implemented');
    elseif A.blockType == 'D' && B.blockType == 'H'
        error('Not been implemented');
    elseif A.blockType == 'H' && B.blockType == 'H'
        assert(size(C.childHMat,1) == size(A.childHMat,1) ...
            && size(C.childHMat,2) == size(B.childHMat,2) ...
            && size(A.childHMat,2) == size(B.childHMat,1) );
        for iti = 1:size(C.childHMat,1)
            for itj = 1:size(C.childHMat,2)
                for itk = 1:size(A.childHMat,2)
                    C.childHMat{iti,itj} = hmtimesh( ...
                        A.childHMat{iti,itk}, B.childHMat{itk,itj}, ...
                        C.childHMat{iti,itj});
                end
            end
        end
    end
end

end

function C = hmtimesd(A,B)

if A.blockType == 'L'
    C = LR2D(A.LRMat*B);
elseif A.blockType == 'D'
    C = A.DMat*B;
elseif A.blockType == 'H'
    hoffset = 0;
    C = zeros(A.height,size(B,2));
    for iti = 1:size(A.childHMat,1)
        woffset = 0;
        for itj = 1:size(A.childHMat,2)
            AC = A.childHMat{iti,itj};
            C(hoffset+(1:AC.height),:) = C(hoffset+(1:AC.height),:) ...
                + hmtimesd(AC,B(woffset+(1:AC.width),:));
            woffset = woffset + AC.width;
        end
        hoffset = hoffset + AC.height;
    end
elseif A.blockType == 'E'
    C = zeros(A.height,size(B,2));
end

end

function C = dmtimesh(A,B)

if B.blockType == 'L'
    C = LR2D(A*B.LRMat);
elseif B.blockType == 'D'
    C = A*B.DMat;
elseif B.blockType == 'H'
    hoffset = 0;
    C = zeros(size(A,1),B.width);
    for itj = 1:size(B.childHMat,2)
        woffset = 0;
        for iti = 1:size(B.childHMat,1)
            BC = B.childHMat{iti,itj};
            C(:,hoffset+(1:BC.width)) = C(:,hoffset+(1:BC.width)) ...
                + dmtimesh(A(:,woffset+(1:BC.height)),BC);
            woffset = woffset + BC.height;
        end
        hoffset = hoffset + BC.width;
    end
elseif B.blockType == 'E'
    C = zeros(size(A,1),B.width);
end

end
