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
MaxRank = 4
minn = 4

Z2C = Z2Cmapper([n,n],minn,reshape(1:n^2,n,n))
A = full(A[Z2C,Z2C]);
@time AH = HMat2dd2h(A, [n,n], [n,n], "WEAK", [0,0], [0,0], 1, EPS, MaxRank, minn)

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

AD = hh2d(AH)
@printf("hh2d error: %.3e\n",norm(AD-A)/norm(A))

BH = hcopy(AH)
CH = hadd(AH,BH)
CD = hh2d(CH)
@printf("hadd error: %.3e\n",norm(CD-A-A)/norm(2*A))

UMat = rand(n^2,MaxRank)
VMat = rand(n^2,MaxRank)
CH = hadd(AH,UMat,VMat)
CD = hh2d(CH)
@printf("hadd error: %.3e\n",norm(CD-A-UMat*VMat')/norm(A+UMat*VMat'))

alpha = 3.1
CH = hscale(alpha,AH)
CD = hh2d(CH)
@printf("hscale error: %.3e\n",norm(CD-alpha*A)/norm(alpha*A))

BH = hcopy(AH)
CH = hmul(AH,BH)
CD = hh2d(CH)
@printf("hmul error: %.3e\n",norm(CD-A*A)/norm(A*A))

Anormfro = hnorm(AH,"fro")
@printf("hnorm fro error: %.3e\n",norm(Anormfro-normfro(A))/normfro(A))
