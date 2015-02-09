include("../src/HMat2d.jl")

function laplacian2d(nx,ny)
    M = speye(nx*ny)
    for x=1:nx
        for y=1:ny
            s = x+(y-1)*nx
            M[s,s] = 4
            if x > 1
                M[s,s-1] = -1
            end
            if x < nx
                M[s,s+1] = -1
            end
            if y > 1
                M[s,s-nx] = -1
            end
            if y < ny
                M[s,s+nx] = -1
            end
        end
    end
    return M
end

n = 8
A = laplacian2d(n,n)

EPS = 1e-2
MaxRank = 2
minn = 4

Z2C = Z2Cmapper([n,n],minn,reshape(1:n^2,n,n))
A = A[Z2C,Z2C];
@time AH = HMat2dd2h(A, [n,n], [n,n], "STANDARD", minn, [0,0], [0,0], 1, EPS, MaxRank, minn)

v = rand(n^2,5)
uext  = A*v
uhmat = hmatvec(AH,v)
@printf("hmatvec error: %.3e\n",norm(uext-uhmat)/norm(uext))

uext  = A'*v
uhmat = hmatTvec(AH,v)
@printf("hmatTvec error: %.3e\n",norm(uext-uhmat)/norm(uext))

uext  = v'*A
uhmat = hvecmat(v',AH)
@printf("hvecmat error: %.3e\n",norm(uext-uhmat)/norm(uext))

uext  = v'*A'
uhmat = hvecmatT(v',AH)
@printf("hvecmatT error: %.3e\n",norm(uext-uhmat)/norm(uext))
