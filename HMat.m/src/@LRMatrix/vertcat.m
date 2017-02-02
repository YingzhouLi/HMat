function A = vertcat( varargin )

height = 0;
rank = 0;
width = varargin{1}.width;
for it = 1:nargin
    height = height + varargin{it}.height;
    rank = rank + size(varargin{it}.UMat,2);
end

UMat = zeros(height,rank);
VMat = zeros(rank,width);

roffset = 0;
hoffset = 0;
for it = 1:nargin
    rank = size(varargin{it}.UMat,2);
    height = varargin{it}.height;
    UMat(hoffset+(1:height),roffset+(1:rank)) = varargin{it}.UMat;
    VMat(:,roffset+(1:rank)) = varargin{it}.VMat;
    roffset = roffset + rank;
    hoffset = hoffset + height;
end

A = LRMatrix(UMat,VMat,varargin{1}.EPS,varargin{1}.MAXRANK);

end