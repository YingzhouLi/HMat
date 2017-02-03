function A = vertcat( varargin )

height = 0;
width = varargin{1}.width;
for it = 1:nargin
    height = height + varargin{it}.height;
end

A = HMatrix;
A.height = height;
A.width = width;
A.blockType = 'H';
A.childHMat = cell(nargin,1);
for it = 1:nargin
    A.childHMat{it} = varargin{it};
end

end