function D = LR2D(A)

if isa(A,'LRMatrix')
    D = A.UMat*A.VMat';
else
    D = A;
end

end