function M = laplacian2d(nx,ny)
    M = speye(nx*ny);
    for x=1:nx
        for y=1:ny
            s = x+(y-1)*nx;
            M(s,s) = 4;
            if x > 1
                M(s,s-1) = -1;
            end
            if x < nx
                M(s,s+1) = -1;
            end
            if y > 1
                M(s,s-nx) = -1;
            end
            if y < ny
                M(s,s+nx) = -1;
            end
        end
    end
end
