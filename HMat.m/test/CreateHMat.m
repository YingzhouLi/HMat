function [varargout]=CreateHMat(varargin)
%getting number of matrices we are going to convert
InputsSize=nargin;
%all the same size, therefore initiating sizes for just 1
[EPS,MaxRank,minn,n1,n2,n3,n4] = InitDefaultHmatVars( varargin{1} );
for i=1:InputsSize
    %creating hmatricies from inputs
    varargin{i}= HMatrix(varargin{i}, [n1,n2],[n3,n4], 'S', [0,0], [0,0],  EPS, MaxRank, minn);
end
%pushing these out
varargout=varargin;
end