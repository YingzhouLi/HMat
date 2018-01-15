addpath('../src');

n = 50;
A = full(laplacian2d(n,n));

EPS = 1e-6;
MaxRank = 2;
minn = 4;

Z2C = Z2Cmapper([n,n],minn,reshape(1:n^2,n,n));
A = full(A(Z2C,Z2C));

tic;
AH = hd2h(A, [n,n], [n,n], 'S', [0,0], [0,0], 1, EPS, MaxRank, minn);
toc;

v = rand(n^2,5);
uext  = A*v;
tic;
uhmat = hmatvec(AH,v);
toc;
fprintf('hmatvec error: %.3e\n',norm(uext-uhmat)/norm(uext));

uext  = A'*v;
tic;
uhmat = hmatTvec(AH,v);
toc;
fprintf('hmatTvec error: %.3e\n',norm(uext-uhmat)/norm(uext));

uext  = v'*A;
tic;
uhmat = hvecmat(v',AH);
toc;
fprintf('hvecmat error: %.3e\n',norm(uext-uhmat)/norm(uext));

uext  = v'*A';
tic;
uhmat = hvecmatT(v',AH);
toc;
fprintf('hvecmatT error: %.3e\n',norm(uext-uhmat)/norm(uext));

AD = hh2d(AH);
fprintf('hh2d error: %.3e\n',norm(AD-A)/norm(A));

BH = hcopy(AH);
CH = haddh(AH,BH);
CD = hh2d(CH);
fprintf('hadd error: %.3e\n',norm(CD-A-A)/norm(2*A));

UMat = rand(n^2,MaxRank);
VMat = rand(n^2,MaxRank);
CH = haddl(AH,UMat,VMat);
CD = hh2d(CH);
fprintf('hadd error: %.3e\n',norm(CD-A-UMat*VMat')/norm(A+UMat*VMat'));

alpha = 3.1;
CH = hscale(alpha,AH);
CD = hh2d(CH);
fprintf('hscale error: %.3e\n',norm(CD-alpha*A)/norm(alpha*A));

BH = hcopy(AH);
tic;
CH = hmul(AH,BH);
toc;
CD = hh2d(CH);
fprintf('hmul error: %.3e\n',norm(CD-A*A)/norm(A*A));

Anormfro = hnormfro(AH);
fprintf('hnorm fro error: %.3e\n',norm(Anormfro-norm(A,'fro'))/norm(A,'fro'));
