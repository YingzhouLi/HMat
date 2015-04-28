function u = hmatTvec(A,v)
    assert(A.height == size(v,1));
    if A.blockType == 'L'
        u = A.VMat*(A.UMat'*v);
    elseif A.blockType == 'D'
        u = A.DMat'*v;
    else
        uoffset = 0;
        u = zeros(A.width,size(v,2));
        for j = 1:4
            voffset = 0;
            for i = 1:4
                C = A.childHMat{i,j};
                u(uoffset+(1:C.width),:) = u(uoffset+(1:C.width),:)...
                        + hmatTvec(C,v(voffset+(1:C.height),:));
                voffset = voffset + C.height;
                if i == 4
                    uoffset = uoffset + C.width;
                end
            end
        end
    end
end

