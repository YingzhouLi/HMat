function u = hmatvec(A,v)
    assert(A.width == size(v,1));
    if A.blockType == 'L'
        u = A.UMat*(A.VMat'*v);
    elseif A.blockType == 'D'
        u = A.DMat*v;
    else
        uoffset = 0;
        u = zeros(A.height,size(v,2));
        for i = 1:4
            voffset = 0;
            for j = 1:4
                C = A.childHMat{i,j};
                u(uoffset+(1:C.height),:) = u(uoffset+(1:C.height),:)...
                        + hmatvec(C,v(voffset+(1:C.width),:));
                voffset = voffset + C.width;
                if j == 4
                    uoffset = uoffset + C.height;
                end
            end
        end
    end
end
