function DH = hd2h(D, nTrg, nSrc, type_admiss, idxTrg, idxSrc, level, EPS, MaxRank, minn)
    height = prod(nTrg);
    width = prod(nSrc);

    if admiss(idxTrg,idxSrc,type_admiss)
        [Utmp,Stmp,Vtmp] = svdtrunc(D,EPS,MaxRank);
        blockType = 'L';
        if max(size(Stmp)) > 0
            UMat = Utmp*sqrt(Stmp);
            VMat = Vtmp*sqrt(Stmp);
        else
            UMat = Utmp;
            VMat = Vtmp;
        end
        DH = struct('height',height,'width',width,...
                    'level',level,...
                    'trg',idxTrg, 'src',idxSrc,...
                    'blockType',blockType,...
                    'UMat',UMat, 'VMat',VMat,...
                    'EPS',EPS, 'MAXRANK',MaxRank, 'MINN',minn);
    elseif max(nTrg) <= minn || max(nSrc) <= minn
        blockType = 'D';
        DH = struct('height',height,'width',width,...
                    'level',level,...
                    'trg',idxTrg, 'src',idxSrc,...
                    'blockType',blockType,...
                    'DMat',full(D),...
                    'EPS',EPS, 'MAXRANK',MaxRank, 'MINN',minn);
    else
        blockType = 'H';
        childHMat = cell(4,4);
        toffset = 0;
        for tx = 0:1
        for ty = 0:1
            trg = 2*idxTrg + [tx,ty];
            txlen = floor(nTrg(1)/2)*(1-tx) + ceil(nTrg(1)/2)*tx;
            tylen = floor(nTrg(2)/2)*(1-ty) + ceil(nTrg(2)/2)*ty;
            tRange = toffset + (1:txlen*tylen);
            toffset = toffset + txlen*tylen;

            soffset = 0;
            for sx = 0:1
            for sy = 0:1
                src = 2*idxSrc + [sx,sy];
                sxlen = floor(nSrc(1)/2)*(1-sx) + ceil(nSrc(1)/2)*sx;
                sylen = floor(nSrc(2)/2)*(1-sy) + ceil(nSrc(2)/2)*sy;
                sRange = soffset + (1:sxlen*sylen);
                soffset = soffset + sxlen*sylen;
                childHMat{tx*2+ty+1,sx*2+sy+1} = hd2h(D(tRange,sRange),...
                            [txlen,tylen],[sxlen,sylen],type_admiss,trg,src,...
                            level+1,EPS,MaxRank,minn);

            end
            end
        end
        end
        DH = struct('height',height,'width',width,...
                    'level',level,...
                    'trg',idxTrg, 'src',idxSrc,...
                    'blockType',blockType,...
                    'childHMat',{childHMat},...
                    'EPS',EPS, 'MAXRANK',MaxRank, 'MINN',minn);
    end
end
