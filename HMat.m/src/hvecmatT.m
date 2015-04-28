function u = hvecmatT(v,A)
    assert(A.width == size(v,2));
    if A.blockType == 'L'
        u = v*A.VMat*A.UMat';
    elseif A.blockType == 'D'
        u = v*A.DMat';
    else
        uoffset = 0;
        u = zeros(size(v,1),A.height);
        for i = 1:4
            voffset = 0;
            for j = 1:4
                C = A.childHMat{i,j};
                u(:,uoffset+(1:C.height)) = u(:,uoffset+(1:C.height))...
                    + hvecmatT(v(:,voffset+(1:C.width)),C);
                voffset = voffset + C.width;
                if j == 4
                    uoffset = uoffset + C.height;
                end
            end
        end
    end
end
