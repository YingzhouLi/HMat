function hNewton(AA::HMat2d,A0::HMat2d,numIter=30)
    tol = 1e-2;
    X = hcopy(A0);
    Xold = hscale(-1,X);
    A = hcopy(AA);
    hscale!(-1,A);
    for iter = 1:numIter
        Z = hmul(X,A);
        hadddiag!(Z,2);
        X = hmul(Z,X);
        EMat = hadd(Xold,X);
        err = hnorm(EMat);
        println(iter,": ",err);
        if err < tol
            break;
        end
        Xold = hscale(-1,X);
    end
    return X;
end
