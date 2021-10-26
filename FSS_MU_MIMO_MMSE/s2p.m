function [y, nz] = s2p(x, np)
%**************************************************************
%S2P    Serial to parallel converter.
%    Y = S2P(X, NP) converts serial data X into NP parallel data.
%    If the length of x is not a multiple of np, zeros are added to fit.
%
%    [Y, NZ] = S2P(X, NP) outputs paralell data and NZ.
%
%   INPUTS
%   =======
%  X= Input data sequence
%  NP= The number of parallel bits
%
%  OUTPUTS
%  =========
%  Y=Output datasequence
%  NZ= The number of zeros added
%
%  Copyright (c) Chang-Jun Ahn 2000 (jun@sasase.ics.keio.ac.jp)
%**************************************************************

[r, c] = size(x);    % x is a r-by-c matrix
if r * c == 0            % if r*c=0, y is empty
    y = [];
    return;
end;
if r == 1                  % if r=1 (i.e., x is a row vector),
    x = x(:);             % convert x into a column vector
    len_x = c;          % length of x = c
else
    len_x = r;           % otherwise length of x = r
end;

nz = mod(len_x, np);
if nz~=0                % If the length of x is not a multiple of np,
    nz=np-nz;
    xadd=zeros(nz,1);
    x=[x;xadd];     % zeros are added to fit
    len_x=len_x+nz;
end;

y=reshape(x, np, len_x/np).';    % Serial to parallel
% y=y(:,np:-1:1);                % reorder from MSB to LSB

