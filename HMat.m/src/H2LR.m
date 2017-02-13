function LRMat = H2LR(A)

warning('H2LR is used for dense matrix');
[Utmp,Stmp,Vtmp] = svdtrunc(A,1e-12,min(size(A)));
LRMat = LRMatrix(Utmp*sqrt(Stmp),Vtmp*sqrt(Stmp),1e-12,min(size(A)));

end