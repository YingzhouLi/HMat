function u = hvecmat(v,A)
    assert(A.height == size(v,2));
    if A.blockType == 'L'
        u = v*A.UMat*A.VMat';
    elseif A.blockType == 'D'
        u = v*A.DMat;
    else
        uoffset = 0;
        u = zeros(size(v,1),A.width);
        for j = 1:4
            voffset = 0;
            for i = 1:4
                C = A.childHMat{i,j};
                u(:,uoffset+(1:C.width)) = u(:,uoffset+(1:C.width))...
                        + hvecmat(v(:,voffset+(1:C.height)),C);
                voffset = voffset + C.height;
                if i == 4
                    uoffset = uoffset + C.width;
                end
            end
        end
    end
end

