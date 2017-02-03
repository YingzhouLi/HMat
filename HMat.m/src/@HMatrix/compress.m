function A = compress(A, mul)

if nargin == 1
    mul = 5;
end

if A.blockType == 'L'
    A.LRMat = compress(A.LRMat,mul);
elseif A.blockType == 'H'
    for i = 1:numel(A.childHMat)
        A.childHMat{i} = compress(A.childHMat{i}, mul);
    end
end

end