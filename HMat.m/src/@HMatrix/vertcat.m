function A = vertcat( varargin )

height = 0;
width = varargin{1}.width;
merge = 0;
hblock = 0;
for it = 1:nargin
    height = height + varargin{it}.height;
    if varargin{it}.blockType == 'H'
        if size(varargin{it}.childHMat,2) == merge || merge == 0
            merge = size(varargin{it}.childHMat,2);
            hblock = hblock + size(varargin{it}.childHMat,1);
        else
            merge = -1;
        end
    else
        merge = -1;
    end
end

if merge == -1
    A = HMatrix;
    A.height = height;
    A.width = width;
    A.blockType = 'H';
    A.childHMat = cell(nargin,1);
    for it = 1:nargin
        A.childHMat{it} = varargin{it};
    end
else
    A = HMatrix;
    A.height = height;
    A.width = width;
    A.blockType = 'H';
    A.childHMat = cell(hblock,merge);
    offset = 0;
    for it = 1:nargin
        for iti = 1:size(varargin{it}.childHMat,1)
            for itj = 1:merge
                A.childHMat{offset + iti,itj} = ...
                    varargin{it}.childHMat{iti,itj};
            end
        end
        offset = offset + size(varargin{it}.childHMat,1);
    end
end

end