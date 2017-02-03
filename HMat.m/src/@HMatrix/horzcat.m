function A = horzcat( varargin )

height = varargin{1}.height;
width = 0;
for it = 1:nargin
    width = width + varargin{it}.width;
end

A = HMatrix;
A.height = height;
A.width = width;
A.blockType = 'H';
A.childHMat = cell(1,nargin);
for it = 1:nargin
    A.childHMat{it} = varargin{it};
end

end