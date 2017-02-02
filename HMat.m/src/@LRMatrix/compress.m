function LR = compress(LR, mul)

if nargin == 1
    mul = 1;
end

if size(LR.UMat,2) > mul*LR.MAXRANK
    [QU,RU] = qr(LR.UMat,0);
    [QV,RV] = qr(LR.VMat,0);
    [Utmp,Stmp,Vtmp] = svdtrunc(RU*RV',LR.EPS,LR.MAXRANK);
    if max(size(Stmp)) > 0
        LR.UMat = QU*(Utmp*sqrt(Stmp));
        LR.VMat = QV*(Vtmp*sqrt(Stmp));
    else
        LR.UMat = QU*Utmp;
        LR.VMat = QV*Vtmp;
    end
end

end