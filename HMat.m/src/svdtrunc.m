function [U,S,V] = svdtrunc(A,eps,mR)
    [U,S,V] = svd(full(A),0);
    if min(size(A))>1 && S(1) > 1e-15
        idx = find(find(diag(S)/S(1)>eps)<=mR);
        U = U(:,idx);
        S = S(idx,idx);
        V = V(:,idx);
    else
        U = zeros(size(A,1),0);
        S = zeros(0);
        V = zeros(size(A,2),0);
    end
end
