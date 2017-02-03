function C = times( A, B )

if isa(A,'HMatrix') && isa(B,'HMatrix')
    C = htimesh(A,B);
elseif isa(A,'HMatrix') && isa(B,'LRMatrix')
    C = htimesl(A,B);
elseif isa(A,'LRMatrix') && isa(B,'HMatrix')
    C = ltimesh(A,B);
else
    C = H2D(A) .* H2D(B);
end

end

function C = htimesh(A,B)

if A.blockType == 'L'
    LRMat = H2L(B);
    C = A;
    C.LRMat = C.LRMat .* LRMat;
    C = compress(C);
elseif A.blockType == 'D'
    C = A;
    C.DMat = C.DMat .* H2D(B);
elseif A.blockType == 'H'
    if B.blockType == 'L'
        C = A .* B.LRMat;
    elseif B.blockType == 'D'
        C = B;
        C.DMat = H2D(A) .* B.DMat;
    else
        assert(all(size(A.childHMat) == size(B.childHMat)));
        C = A;
        for iti = 1:size(A.childHMat,1)
            for itj = 1:size(A.childHMat,2)
                C.childHMat{iti,itj} = htimesh( A.childHMat{iti,itj}, ...
                    B.childHMat{iti,itj} );
            end
        end
    end
end

end

function C = htimesl(A,LRMat)

C = A;
if A.blockType == 'L'
    C.LRMat = A.LRMat .* LRMat;
    C = compress(C);
elseif A.blockType == 'D'
    C.DMat = A.DMat .* LR2D(LRMat);
elseif A.blockType == 'H'
    hoffset = 0;
    for iti = 1:size(A.childHMat,1)
        woffset = 0;
        for itj = 1:size(A.childHMat,2)
            height = A.childHMat{iti,itj}.height;
            width = A.childHMat{iti,itj}.width;
            C.childHMat{iti,itj} = htimesl( A.childHMat{iti,itj},...
                LRMat(hoffset+(1:height), woffset+(1:width)));
            woffset = woffset + width;
        end
        hoffset = hoffset + height;
    end
end

end

function C = ltimesh(LRMat,A)

C = A;
if A.blockType == 'L'
    C.LRMat = A.LRMat .* LRMat;
    C = compress(C);
elseif A.blockType == 'D'
    C.DMat = A.DMat .* LR2D(LRMat);
elseif A.blockType == 'H'
    hoffset = 0;
    for iti = 1:size(A.childHMat,1)
        woffset = 0;
        for itj = 1:size(A.childHMat,2)
            height = A.childHMat{iti,itj}.height;
            width = A.childHMat{iti,itj}.width;
            C.childHMat{iti,itj} = ltimesh( ...
                LRMat(hoffset+(1:height), woffset+(1:width)), ...
                A.childHMat{iti,itj});
            woffset = woffset + width;
        end
        hoffset = hoffset + height;
    end
end

end