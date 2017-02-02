function C = hempty(A)

C = A;
if A.blockType == 'L'
    C.UMat = zeros(A.height,0);
    C.VMat = zeros(A.width,0);
elseif A.blockType == 'D'
    C.DMat = zeros(A.height,A.width);
elseif A.blockType == 'H'
    C.childHMat = cell(size(A.childHMat));
    for i = 1:length(A.childHMat)
        C.childHMat{i} = hempty(A.childHMat{i});
    end
end

end