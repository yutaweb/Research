function d = b2d(b, p)
%*******************************************************************
%B2D    Converts binary to decimal format.
%  D = B2D(B) converts a binary representation vector B to a decimal
%  presentation D. When B is a matrix, the output D is a column vector
%  with each element being the transfer of a row of B.
%  The first element in B represents the highest binary bit. For example
%  b2d([1 0]) results in 2; b2d([0 1]) results in 1.
%
%  D = B2D(B, P) converts a P-based vector to a decimal.
%
%  Copyright (c) Chang-Jun Ahn 2000 (jun@sasase.ics.keio.ac.jp) 
%********************************************************************

%********************************************************************
% check the special cases.
%********************************************************************
[n,m] = size(b);
if min([m,n]) < 1
    d = [];
    return;
elseif min([n,m]) == 1
    b = b(:)';
    m = max([n,m]);
    n = 1;
end;

if nargin < 2
    p = 2;
end;
d = b * p.^[m-1:-1:0]';