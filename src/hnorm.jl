function hnorm(A::HMat2d,p=Inf)
    if p == Inf
        return hnorminf(A)
    elseif p == 1
        error("Not been implemented")
    elseif p == 2
        error("Not been implemented")
    else
        error("hnorm only support p = 1,2,Inf")
    end
end