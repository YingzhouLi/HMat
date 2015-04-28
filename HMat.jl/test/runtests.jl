include("../src/HMat2d.jl");

function laplacian2d(nx,ny)
    M = speye(nx*ny);
    for x=1:nx
        for y=1:ny
            s = x+(y-1)*nx;
            M[s,s] = 4;
            if x > 1
                M[s,s-1] = -1;
            end
            if x < nx
                M[s,s+1] = -1;
            end
            if y > 1
                M[s,s-nx] = -1;
            end
            if y < ny
                M[s,s+nx] = -1;
            end
        end
    end
    return M;
end

n = 32;
A = full(laplacian2d(n,n));

EPS = 1e-6;
MaxRank = 8;
minn = 4;

Z2C = Z2Cmapper([n,n],minn,reshape(1:n^2,n,n));
A = full(A[Z2C,Z2C]);

@time AH = HMat2dd2h(A, [n,n], [n,n], STANDARD, [0,0], [0,0], 1, EPS, MaxRank, minn);

v = rand(n^2,5);
uext  = A*v;
@time uhmat = hmatvec(AH,v);
@time uhmat = hmatvec(AH,v);
@printf("hmatvec error: %.3e\n",norm(uext-uhmat)/norm(uext));

uext  = A'*v;
@time uhmat = hmatTvec(AH,v);
@time uhmat = hmatTvec(AH,v);
@printf("hmatTvec error: %.3e\n",norm(uext-uhmat)/norm(uext));

uext  = v'*A;
@time uhmat = hvecmat(v',AH);
@time uhmat = hvecmat(v',AH);
@printf("hvecmat error: %.3e\n",norm(uext-uhmat)/norm(uext));

uext  = v'*A';
@time uhmat = hvecmatT(v',AH);
@time uhmat = hvecmatT(v',AH);
@printf("hvecmatT error: %.3e\n",norm(uext-uhmat)/norm(uext));

@time AD = hh2d(AH);
@printf("hh2d error: %.3e\n",norm(AD-A)/norm(A));

BH = hcopy(AH);
@time CH = hadd(AH,BH);
@time CH = hadd(AH,BH);
CD = hh2d(CH);
@printf("hadd error: %.3e\n",norm(CD-A-A)/norm(2*A));

UMat = rand(n^2,MaxRank);
VMat = rand(n^2,MaxRank);
@time CH = hadd(AH,UMat,VMat);
CD = hh2d(CH);
@printf("hadd error: %.3e\n",norm(CD-A-UMat*VMat')/norm(A+UMat*VMat'));

alpha = 3.1;
@time CH = hscale(alpha,AH);
CD = hh2d(CH);
@printf("hscale error: %.3e\n",norm(CD-alpha*A)/norm(alpha*A));

BH = hcopy(AH);
@time CH = hmul(AH,BH);
CD = hh2d(CH);
AA = A*A;
@printf("hmul error: %.3e\n",norm(CD-A*A)/norm(A*A));
Err = CD-A*A;

@time Anormfro = hnorm(AH,"fro");
@printf("hnorm fro error: %.3e\n",norm(Anormfro-vecnorm(A))/vecnorm(A));
