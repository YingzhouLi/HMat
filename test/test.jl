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

EPS = 1e-6
MaxRank = 2
minn = 8

for n = [16 32 64]
    A = full(laplacian2d(n,n))
    Z2C = Z2Cmapper([n,n],minn,reshape(1:n^2,n,n))
    A = full(A[Z2C,Z2C]);
    AH = HMat2dd2h(A, [n,n], [n,n], WEAK, [0,0], [0,0], 1, EPS, MaxRank, minn)
    BH = hcopy(AH)
    v = rand(n^2,5)
    tic()
    uhmat = hmatvec(AH,v)
    toc()
    tic()
    CH = hmul(AH,BH)
    toc()
    CD = hh2d(CH)
    @printf("N^2 = %4d^2, hmul error: %.3e\n",n,norm(CD-A*A)/norm(A*A))
end
