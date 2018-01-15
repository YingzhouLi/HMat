export LRMat;

struct LRMat{T<:Number}
    # variables
    height::    Int
    width::     Int
    UMat::      Array{T,2}
    VMat::      Array{T,2}
    # global settings
    EPS::       Float64
    MAXRANK::   Int

    function LRMat(D,Eps,MaxRank)
        h = size(D,1);
        w = size(D,2);
        [U,S,V] = svdtrunc(D,Eps,MaxRank);
        new(h,w,U*sqrt(S),V*sqrt(S),Eps,MaxRank);
    end
end
