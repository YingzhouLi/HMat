function A = compress(A, mul)

if nargin == 1
    mul = 5;
end

if A.blockType == 'L'
    if size(A.UMat,2) > mul*A.MAXRANK
        [QU,RU] = qr(A.UMat,0);
        [QV,RV] = qr(A.VMat,0);
        [Utmp,Stmp,Vtmp] = svdtrunc(RU*RV',A.EPS,A.MAXRANK);
        if max(size(Stmp)) > 0
            A.UMat = QU*(Utmp*sqrt(Stmp));
            A.VMat = QV*(Vtmp*sqrt(Stmp));
        else
            A.UMat = QU*Utmp;
            A.VMat = QV*Vtmp;
        end
    end
elseif A.blockType == 'H'
    for i = 1:length(A.childHMat)
        A.childHMat{i} = compress(A.childHMat{i}, mul);
    end
end

end