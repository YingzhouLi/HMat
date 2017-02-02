classdef LRMatrix
    properties
        height = 0
        width = 0
        UMat
        VMat
        EPS
        MAXRANK
    end
    
    methods
        
        function LR = LRMatrix( varargin )
            
            if nargin == 3
                
                D       = varargin{1};
                EPS     = varargin{2};
                MaxRank = varargin{3};
                
                [LR.height,LR.width] = size(D);
                LR.EPS = EPS;
                LR.MAXRANK = MaxRank;
                
                [Utmp,Stmp,Vtmp] = svdtrunc(D,EPS,MaxRank);
                if max(size(Stmp)) > 0
                    LR.UMat = Utmp*sqrt(Stmp);
                    LR.VMat = Vtmp*sqrt(Stmp);
                else
                    LR.UMat = Utmp;
                    LR.VMat = Vtmp;
                end
                
            elseif nargin == 4
                
                LR.UMat    = varargin{1};
                LR.VMat    = varargin{2};
                LR.EPS     = varargin{3};
                LR.MAXRANK = varargin{4};
                
                LR.height  = size(LR.UMat,1);
                LR.width   = size(LR.VMat,1);
                
            end
            
        end
        
        LR = compress(LR, mul);
        
    end
end