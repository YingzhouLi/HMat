function A = horzcat( varargin )

height = varargin{1}.height;
width = 0;
merge = 0;
wblock = 0;
for it = 1:nargin
    width = width + varargin{it}.width;
    if varargin{it}.blockType == 'H'
        if size(varargin{it}.childHMat,1) == merge || merge == 0
            merge = size(varargin{it}.childHMat,1);
            wblock = wblock + size(varargin{it}.childHMat,2);
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
A.childHMat = cell(1,nargin);
for it = 1:nargin
    A.childHMat{it} = varargin{it};
end
else
    A = HMatrix;
    A.height = height;
    A.width = width;
    A.blockType = 'H';
    A.childHMat = cell(merge,wblock);
    offset = 0;
    for it = 1:nargin
        for iti = 1:merge
            for itj = 1:size(varargin{it}.childHMat,2)
                A.childHMat{iti,offset + itj} = ...
                    varargin{it}.childHMat{iti,itj};
            end
        end
        offset = offset + size(varargin{it}.childHMat,2);
    end
end

end