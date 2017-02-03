addpath('../src');

EPS = 1e-6;
MaxRank = 2;
minn = 4;

for n = 32
    A = full(laplacian2d(n,n));

    Z2C = Z2Cmapper([n,n],minn,reshape(1:n^2,n,n));
    A = full(A(Z2C,Z2C));

    tic;
    AH = HMatrix(A, [n,n], [n,n], 'S', [0,0], [0,0], EPS, MaxRank, minn);
    toc;

    tic;
	BH = AH;
    toc;
    
    tic;
    CH = AH*BH;
    toc;
    %CD = hh2d(CH);
    %fprintf('N^2 = %3d^2, hmul error: %.3e\n',n,norm(CD-A*A)/norm(A*A));
end
