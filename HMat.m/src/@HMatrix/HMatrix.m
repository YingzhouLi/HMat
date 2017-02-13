classdef HMatrix
    properties
        height = 0
        width = 0
        blockType = 'E'
        LRMat
        DMat
        childHMat
        EPS
        MAXRANK
    end
    
    methods
        
        function DH = HMatrix( D, nTrg, nSrc, type_admiss, idxTrg, idxSrc, ...
                EPS, MaxRank, minn)
            
            if nargin == 0
                return;
            end
            
            DH.height = prod(nTrg);
            DH.width = prod(nSrc);
            DH.EPS = EPS;
            DH.MAXRANK = MaxRank;
            
            if admiss(idxTrg,idxSrc,type_admiss)
                DH.blockType = 'L';
                DH.LRMat = LRMatrix(D,EPS,MaxRank);
            elseif max(nTrg) <= minn || max(nSrc) <= minn
                DH.blockType = 'D';
                DH.DMat = full(D);
            else
                DH.blockType = 'H';
                
                nT = 2*(max(nTrg) > 1);
                nS = 2*(max(nSrc) > 1);
                DH.childHMat = cell(nT,nS);
                
                toffset = 0;
                for itT = 0:nT-1
                    [~,idxT] = max(nTrg);
                    trg = idxTrg;
                    trg(idxT) = idxTrg(idxT)*nT + itT;
                    tlen = nTrg;
                    tlen(idxT) = floor(nTrg(idxT)/nT)*(1-itT) ...
                        + ceil(nTrg(idxT)/nT)*itT;
                    tRange = toffset + (1:prod(tlen));
                    toffset = toffset + prod(tlen);
                    
                    soffset = 0;
                    for itS = 0:nS-1
                        [~,idxS] = max(nSrc);
                        src = idxSrc;
                        src(idxS) = idxSrc(idxS)*nS + itS;
                        slen = nSrc;
                        slen(idxS) = floor(nSrc(idxS)/nS)*(1-itS) ...
                            + ceil(nSrc(idxS)/nS)*itS;
                        sRange = soffset + (1:prod(slen));
                        soffset = soffset + prod(slen);
                        
                        DH.childHMat{itT+1,itS+1} = ...
                            HMatrix( D(tRange, sRange), ...
                            tlen, slen, type_admiss, trg, src, ...
                            EPS, MaxRank, minn );
                        
                    end
                end
            end
        end
        
        A = compress(A, mul);
        
    end
end
