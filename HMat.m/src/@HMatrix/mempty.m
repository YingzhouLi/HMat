function C = mempty(A,B)

C = HMatrix();
C.height = A.height;
C.width = B.width;
C.EPS = A.EPS;
C.MAXRANK = A.MAXRANK;

if A.blockType == 'E' || B.blockType == 'E'
    C.blockType = 'E';
elseif A.blockType == 'D' || B.blockType == 'D'
    C.blockType = 'D';
    C.DMat = zeros(A.height,B.width);
elseif A.blockType == 'L' || B.blockType == 'L'
    C.blockType = 'L';
    C.LRMat = LRMatrix( zeros(A.height,0), zeros(B.width,0), ...
        A.EPS, A.MAXRANK );
elseif A.blockType == 'H' && B.blockType == 'H'
    C.blockType = 'H';
    C.childHMat = cell(size(A.childHMat,1),size(B.childHMat,2));
    for iti = 1:size(C.childHMat,1)
        for itj = 1:size(C.childHMat,2)
            complvl = 0;
            for itk = 1:size(A.childHMat,2)
                T1 = A.childHMat{iti,itk}.blockType;
                T2 = B.childHMat{itk,itj}.blockType;
                if complexlevel(T1,T2) > complvl
                    kk = itk;
                end
            end
            C.childHMat{iti,itj} = ...
                mempty(A.childHMat{iti,kk},B.childHMat{kk,itj});
        end
    end
end

end

function lvl = complexlevel(T1,T2)

lvl = 5;

if T1 == 'E' || T2 == 'E'
    lvl = 1;
    return;
end
if T1 == 'D' || T2 == 'D'
    lvl = 2;
    return;
end
if T1 == 'L' || T2 == 'L'
    lvl = 3;
    return;
end
if T1 == 'H' && T2 == 'H'
    lvl = 4;
    return;
end

end