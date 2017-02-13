function C = empty(A,sig)

if nargin < 2
    sig = 'X';
end

if sig == 'E'
    C = HMatrix();
    C.height = A.height;
    C.width = A.width;
    C.blockType = 'E';
    C.EPS = A.EPS;
    C.MAXRANK = A.MAXRANK;
    return;
end

if A.blockType == 'L'
    C = A;
    C.LRMat = LRMatrix( zeros(A.height,0), zeros(A.width,0), ...
        A.LRMat.EPS, A.LRMat.MAXRANK );
elseif A.blockType == 'D'
    C = A;
    C.DMat = zeros(A.height,A.width);
elseif A.blockType == 'H'
    if sig == 'L' || sig == 'U'
        C = HMatrix();
        C.height = A.height;
        C.width = A.width;
        C.blockType = 'H';
        C.childHMat = cell(size(A.childHMat));
        C.EPS = A.EPS;
        C.MAXRANK = A.MAXRANK;
        for iti = 1:size(A.childHMat,1)
            for itj = 1:size(A.childHMat,2)
                if iti == itj
                    C.childHMat{iti,itj} = empty(A.childHMat{iti,itj},sig);
                elseif ( iti < itj && sig == 'U' ) ...
                        || ( iti > itj && sig == 'L' )
                    C.childHMat{iti,itj} = empty(A.childHMat{iti,itj});
                else
                    C.childHMat{iti,itj} = empty(A.childHMat{iti,itj},'E');
                end
            end
        end
    else
        C = A;
        for it = 1:numel(A.childHMat)
            C.childHMat{it} = empty(A.childHMat{it});
        end
    end
else
    C = A;
end

end