function A = horzcat( varargin )

height = varargin{1}.height;
rank = 0;
width = 0;
for it = 1:nargin
    width = width + varargin{it}.width;
    rank = rank + size(varargin{it}.UMat,2);
end

UMat = zeros(height,rank);
VMat = zeros(rank,width);

roffset = 0;
woffset = 0;
for it = 1:nargin
    rank = size(varargin{it}.UMat,2);
    width = varargin{it}.width;
    UMat(:,roffset+(1:rank)) = varargin{it}.UMat;
    VMat(woffset+(1:width),roffset+(1:rank)) = varargin{it}.VMat;
    roffset = roffset + rank;
    woffset = woffset + width;
end

A = LRMatrix(UMat,VMat,varargin{1}.EPS,varargin{1}.MAXRANK);

end