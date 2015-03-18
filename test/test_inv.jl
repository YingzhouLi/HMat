include("../src/HMat2d.jl");

function laplacian2d(DA,DV,nx,ny,h)
    M = speye(nx*ny);
    for x=1:nx
        for y=1:ny
            s = x+(y-1)*nx;
            M[s,s] = DV[x,y];
            if x > 1
                M[s,s-1] = -(DA[x,y+1]+DA[x+1,y+1])/2/h^2;
            end
            M[s,s] += (DA[x,y+1]+DA[x+1,y+1])/2/h^2;
            if x < nx
                M[s,s+1] = -(DA[x+2,y+1]+DA[x+1,y+1])/2/h^2;
            end
            M[s,s] += (DA[x+2,y+1]+DA[x+1,y+1])/2/h^2;
            if y > 1
                M[s,s-nx] = -(DA[x+1,y]+DA[x+1,y+1])/2/h^2;
            end
            M[s,s] += (DA[x+1,y]+DA[x+1,y+1])/2/h^2;
            if y < ny
                M[s,s+nx] = -(DA[x+1,y+2]+DA[x+1,y+1])/2/h^2;
            end
            M[s,s] += (DA[x+1,y+2]+DA[x+1,y+1])/2/h^2;
        end
    end
    return M;
end

EPS = 1e-6;
MaxRank = 4;
minn = 4;

for iter = 1:10

    n=32;
    h=1/n;
    DA = abs(randn(n+2,n+2))+0.01;
    DV = zeros(n,n);
    A = full(laplacian2d(DA,DV,n,n,h));
    Z2C = Z2Cmapper([n,n],minn,reshape(1:n^2,n,n));
    A = inv(full(A[Z2C,Z2C]));
    AHS = HMat2dd2h(A, [n,n], [n,n], STANDARD, [0,0], [0,0], 1, EPS, MaxRank, minn);
    AHSinv = hcopy(AHS);
    hinv!(AHSinv);
    ASinv = hh2d(AHSinv);
    errS = ASinv*A-eye(n^2);
    @printf("Direct inversion of Standard HMat error: %.2e\n",norm(errS));

    X0 = hidentity(AHS);
    hscale!(1/norm(A),X0);
    AHNinv = hNewton(AHS,X0,18);
    ANinv = hh2d(AHNinv);
    errN = ANinv*A-eye(n^2);
    @printf("Newton inversion of Standard HMat error1: %.2e\n\n",norm(errN));

    AHW = HMat2dd2h(A, [n,n], [n,n], WEAK, [0,0], [0,0], 1, EPS, MaxRank, minn);
    AHWinv = hcopy(AHW);
    hinv!(AHWinv);
    Dtmp = hh2d(AHWinv);
    X0 = HMat2dd2h(Dtmp, [n,n], [n,n], STANDARD, [0,0], [0,0], 1, EPS, MaxRank, minn);
    AHNinv = hNewton(AHS,X0,10);
    ANinv = hh2d(AHNinv);
    errN = ANinv*A-eye(n^2);
    @printf("Newton inversion of Standard HMat error2: %.2e\n\n",norm(errN));

end
