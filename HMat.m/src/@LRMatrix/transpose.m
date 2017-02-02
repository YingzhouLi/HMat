function A = transpose( A )

height = A.height;
A.height = A.width;
A.width = height;

UMat = A.UMat;
A.UMat = conj(A.VMat);
A.VMat = conj(UMat);

end