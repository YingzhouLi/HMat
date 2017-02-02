function A = ctranspose( A )

height = A.height;
A.height = A.width;
A.width = height;

UMat = A.UMat;
A.UMat = A.VMat;
A.VMat = UMat;

end