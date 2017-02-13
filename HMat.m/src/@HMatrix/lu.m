function [L,U] = lu( A )

L = empty(A,'L');
U = empty(A,'U');

if A.blockType == 'L'
    error('A cannot be low rank matrix');
elseif A.blockType == 'D'
    [L.DMat,U.DMat] = lu(A.DMat);
elseif A.blockType == 'H'
    for iti = 1:min(size(A.childHMat))
        [L.childHMat{iti,iti},U.childHMat{iti,iti}] = ...
            lu(A.childHMat{iti,iti});
        % Factor L
        for itj = iti+1:size(A.childHMat,1)
            L.childHMat{itj,iti} = ...
                A.childHMat{itj,iti}/U.childHMat{iti,iti};
        end
        % Factor U
        for itj = iti+1:size(A.childHMat,2)
            U.childHMat{iti,itj} = ...
                L.childHMat{iti,iti}\A.childHMat{iti,itj};
        end
        % Update A
        for itj = iti+1:size(A.childHMat,1)
            for itk = iti+1:size(A.childHMat,2)
                A.childHMat{itj,itk} = A.childHMat{itj,itk} ...
                    - L.childHMat{itj,iti}*U.childHMat{iti,itk};
            end
        end
    end
end

end