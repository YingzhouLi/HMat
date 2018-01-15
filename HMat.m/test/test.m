addpath('../src');

EPS = 1e-6;
MaxRank = 12;
minn = 4;

for n = 50
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
    CD = H2D(CH);
    fprintf('N^2 = %3d^2, hmul error: %.3e\n',n,norm(CD-A*A)/norm(A*A));
    
    x = randn(size(A,2),1);
    y = A*x;
    
    tic;
    xest = AH\y;
    toc;
    fprintf('H matrix solve error: %.3e\n',norm(x-xest)/norm(x));
    
    tic;
    [LH,UH] = lu(AH);
    toc;
    
    fprintf('lu factorization error: %.3e\n', ...
        norm(H2D(LH)*H2D(UH) - A)/norm(A));
    
    xest = UH\(LH\y);
    fprintf('lu solve error: %.3e\n',norm(x-xest)/norm(x));
    
end