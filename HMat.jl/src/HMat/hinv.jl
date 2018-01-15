function hinv!(A::HMat)
    assert( A.height == A.width );
    if A.blockType == LOWRANK
        error("Mistake in the code");
    elseif A.blockType == DENSE
        A.DMat = inv(A.DMat);
    elseif A.blockType == HMAT
        B = hcopy(A);
        hidentity!(A);
        for l = 1:4
            A.childHMat[l,l] = hcopy(B.childHMat[l,l]);
            hinv!(A.childHMat[l,l]);
            for j = 1:(l-1)
                C = hcopy(A.childHMat[l,j]);
                A.childHMat[l,j] = hmul(A.childHMat[l,l],C);
            end
            for j = (l+1):4
                C = hcopy(B.childHMat[l,j]);
                B.childHMat[l,j] = hmul(A.childHMat[l,l],C);
            end
            for i = (l+1):4
                for j = 1:l
                    C = hscale(-1,A.childHMat[l,j]);
                    hmul(B.childHMat[i,l],C,A.childHMat[i,j]);
                end
                for j = (l+1):4
                    C = hscale(-1,B.childHMat[l,j]);
                    hmul(B.childHMat[i,l],C,B.childHMat[i,j]);
                end
            end
        end
        for l = 4:-1:1
            for i = (l-1):-1:1
                for j = 1:4
                    C = hscale(-1,A.childHMat[l,j]);
                    hmul(B.childHMat[i,l],C,A.childHMat[i,j]);
                end
            end
        end
    end
end

