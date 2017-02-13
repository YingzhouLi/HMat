function C = mldivide( A, B )

if isa(A,'HMatrix') && isa(B,'HMatrix')
    C = mempty(A',B);
    C = hmldivideh(A,B,C);
elseif isa(A,'HMatrix') && isa(B,'LRMatrix')
    C = LRMatrix(A\B.UMat,B.VMat,B.EPS,B.MAXRANK);
elseif isa(A,'HMatrix')
    C = hmldivided(A,B);
elseif isa(B,'HMatrix')
    C = A\H2D(B);
end

end

function C = hmldivideh(A,B,C)

if C.blockType == 'L'
    if A.blockType == 'L'
        error('A cannot be low rank matrix');
    elseif ( A.blockType == 'D' || A.blockType == 'H' ) ...
            && B.blockType == 'L'
        C.LRMat = A\B.LRMat + C.LRMat;
    else
        Rcol = randn(B.width,C.LRMat.MAXRANK+5);
        Rrow = randn(A.height,C.LRMat.MAXRANK+5);
        ABRcol = A\(B*Rcol);
        BTATRrow = B'*(A'\Rrow);
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
    C.DMat = C.DMat + H2D(A)\H2D(B);
else
    if A.blockType == 'L'
        error('A cannot be low rank matrix');
    elseif ( A.blockType == 'D' || A.blockType == 'H' ) ...
            && B.blockType == 'L'
        C = C + A\B.LRMat;
    elseif A.blockType == 'H' && B.blockType == 'D'
        C = C + A\B.DMat;
    elseif A.blockType == 'D' && B.blockType == 'H'
        C = C + A.DMat\B;
    elseif A.blockType == 'H' && B.blockType == 'H'
        if islower(A)
            for iti = 1:size(C.childHMat,1)
                for itj = 1:size(C.childHMat,2)
                    Btmp = B.childHMat{iti,itj};
                    for itk = 1:iti-1
                        Btmp = Btmp ...
                            - A.childHMat{iti,itk}*C.childHMat{itk,itj};
                    end
                    C.childHMat{iti,itj} = hmldivideh( ...
                        A.childHMat{iti,iti}, Btmp, ...
                        C.childHMat{iti,itj});
                end
            end
        elseif isupper(A)
            for iti = size(C.childHMat,1):-1:1
                for itj = size(C.childHMat,2):-1:1
                    Btmp = B.childHMat{iti,itj};
                    for itk = size(C.childHMat,1):-1:iti+1
                        Btmp = Btmp ...
                            - A.childHMat{iti,itk}*C.childHMat{itk,itj};
                    end
                    C.childHMat{iti,itj} = hmldivideh( ...
                        A.childHMat{iti,iti}, Btmp, ...
                        C.childHMat{iti,itj});
                end
            end
        else
            [L,U] = lu(A);
            C = C + U\(L\B);
        end
    end
end

end

function C = hmldivided(A,B)

if A.blockType == 'L'
    error('A cannot be low rank matrix');
elseif A.blockType == 'D'
    C = A.DMat\B;
else
    C = zeros(A.width,size(B,2));
    if islower(A)
        hoffset = 0;
        for iti = 1:size(A.childHMat,1)
            AC = A.childHMat{iti,iti};
            Btmp = B(hoffset+(1:AC.height),:);
            offset = 0;
            for itj = 1:iti-1
                Cidx = offset + (1:A.childHMat{iti,itj}.width);
                Btmp = Btmp ...
                    - A.childHMat{iti,itj}*C(Cidx,:);
                offset = Cidx(end);
            end
            C(hoffset+(1:AC.width),:) = AC \ Btmp;
            hoffset = hoffset + AC.height;
        end
    elseif isupper(A)
        hoffset = A.height+1;
        for iti = size(A.childHMat,1):-1:1
            AC = A.childHMat{iti,iti};
            Btmp = B(hoffset-(AC.height:-1:1),:);
            offset = A.width+1;
            for itj = size(A.childHMat,1):-1:iti+1
                Cidx = offset - (A.childHMat{iti,itj}.width:-1:1);
                Btmp = Btmp ...
                    - A.childHMat{iti,itj}*C(Cidx,:);
                offset = Cidx(1);
            end
            C(hoffset-(AC.width:-1:1),:) = AC \ Btmp;
            hoffset = hoffset - AC.height;
        end
    else
        [L,U] = lu(A);
        C = C + U\(L\B);
    end
end

end
