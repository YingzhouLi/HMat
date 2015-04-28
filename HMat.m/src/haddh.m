function C = haddh(A,B)
    assert(A.height == B.height)
    assert(A.width == B.width)

    if A.blockType == 'L'
        if B.blockType == 'L'
            UMat = [A.UMat B.UMat];
            VMat = [A.VMat B.VMat];
        else
            [Utmp,Vtmp] = hh2l(B);
            UMat = [A.UMat Utmp];
            VMat = [A.VMat Vtmp];
        end
        C = struct('height',A.height,'width',A.width,...
                    'level',A.level,...
                    'trg',A.trg, 'src',A.src,...
                    'blockType',A.blockType,...
                    'UMat',UMat, 'VMat',VMat,...
                    'EPS',A.EPS, 'MAXRANK',A.MAXRANK, 'MINN',A.MINN);
        C = hcompress(C);
    elseif A.blockType == 'D'
        DMat = A.DMat + hh2d(B);
        C = struct('height',A.height,'width',A.width,...
                    'level',A.level,...
                    'trg',A.trg, 'src',A.src,...
                    'blockType',A.blockType,...
                    'DMat',DMat,...
                    'EPS',A.EPS, 'MAXRANK',A.MAXRANK, 'MINN',A.MINN);
    elseif A.blockType == 'H'
        if B.blockType == 'L'
            C = haddl(A,B.UMat,B.VMat);
        else
            childHMat = cell(4,4);
            if B.blockType == 'D'
                hoffset = 0;
                for i = 1:4
                    woffset = 0;
                    for j = 1:4
                        BHheight = A.childHMat{i,j}.height;
                        BHwidth = A.childHMat{i,j}.width;
                        BH = struct('height',BHheight,...
                                    'width',BHwidth,...
                                    'blockType','D',...
                                    'DMat',B.DMat(hoffset+(1:BHheight),...
                                                  woffset+(1:BHwidth)));
                        woffset = woffset + BHwidth;
                        if j == 4
                            hoffset = hoffset + BHheight;
                        end
                        childHMat{i,j} = haddh(A.childHMat{i,j},BH);
                    end
                end
            else
                for i = 1:4
                    for j = 1:4
                        childHMat{i,j} = haddh(A.childHMat{i,j},B.childHMat{i,j});
                    end
                end
            end
            C = struct('height',A.height,'width',A.width,...
                        'level',A.level,...
                        'trg',A.trg, 'src',A.src,...
                        'blockType',A.blockType,...
                        'childHMat',{childHMat},...
                        'EPS',A.EPS, 'MAXRANK',A.MAXRANK, 'MINN',A.MINN);
        end
    end
end
