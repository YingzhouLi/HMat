function LRMat = H2L(A)

Rcol = randn(A.width,A.MAXRANK+5);
Rrow = randn(A.height,A.MAXRANK+5);
ARcol = A*Rcol;
ATRrow = A'*Rrow;
[Qcol,~] = qr(ARcol,0);
[Qrow,~] = qr(ATRrow,0);
[Utmp,Stmp,Vtmp] = svdtrunc(pinv(Rrow'*Qcol)*Rrow'*ARcol*pinv(Qrow'*Rcol),...
    A.EPS,A.MAXRANK);
if max(size(Stmp)) > 0;
    UMat = Qcol*(Utmp*sqrt(Stmp));
    VMat = Qrow*(Vtmp*sqrt(Stmp));
else
    UMat = Qcol*Utmp;
    VMat = Qrow*Vtmp;
end
LRMat = LRMatrix(UMat,VMat,A.EPS,A.MAXRANK);

end