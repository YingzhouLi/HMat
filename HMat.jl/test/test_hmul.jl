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

EPS = 1e-6;
MaxRank = 2;
minn = 4;

n=16;
A = full(laplacian2d(n,n));
Z2C = Z2Cmapper([n,n],minn,reshape(1:n^2,n,n));
A = inv(full(A[Z2C,Z2C]));
AH16 = HMat2dd2h(A, [n,n], [n,n], EDGE, [0,0], [0,0], 1, EPS, MaxRank, minn);

n=32;
A = full(laplacian2d(n,n));
Z2C = Z2Cmapper([n,n],minn,reshape(1:n^2,n,n));
A = inv(full(A[Z2C,Z2C]));
AH32 = HMat2dd2h(A, [n,n], [n,n], EDGE, [0,0], [0,0], 1, EPS, MaxRank, minn);

n=64;
A = full(laplacian2d(n,n));
Z2C = Z2Cmapper([n,n],minn,reshape(1:n^2,n,n));
A = inv(full(A[Z2C,Z2C]));
AH64 = HMat2dd2h(A, [n,n], [n,n], EDGE, [0,0], [0,0], 1, EPS, MaxRank, minn);

n=16;
BH16 = hcopy(AH16);
CH16 = hmul(AH16,BH16);
tic();
for i=1:10
CH16 = hmul(AH16,BH16);
end
toc();

n=32;
BH32 = hcopy(AH32);
CH32 = hmul(AH32,BH32);
tic();
for i=1:10
CH32 = hmul(AH32,BH32);
end
toc();

n=64;
BH64 = hcopy(AH64);
CH64 = hmul(AH64,BH64);
tic();
for i=1:10
CH64 = hmul(AH64,BH64);
end
toc();
