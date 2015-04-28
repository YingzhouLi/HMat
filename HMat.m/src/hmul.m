function C = hmul(A,B)
    assert(A.width == B.height);
    C = hempty(A);
    C = hmuladd(A,B,C);
    C = hcompress(C,1);
end
