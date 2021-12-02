function b = d2b(d, n, p)
%************************************************************************
%D2B    Converts a positive decimal number to a binary number.
%   B = D2B(D) converts a positive integer decimal vector D to a binary
%   matrix B. Each row of the binary matrix B represents the corresponding
%   number in D.
%
%   B = D2B(D, N) specifies the number of column for matrix B. N is a
%   positive scalar number.
%
%   B = D2B(D, N, P) converts a decimal vector D to base P matrix.
%
%   The first element in B represents the lowest binary bit. For example
%   D2B(1, 2) results in [0 1]; D2B(2, 2) results in [1 0].
%
%   Copyright (c) Chang-Jun Ahn 2000 (jun@sasase.ics.keio.ac.jp)
%************************************************************************

%************************************************************************
% routine check
%************************************************************************
d = d(:);
len_d = length(d);
if min(d) < 0
    error('Cannot convert a negative number');
elseif ~isempty(find(d==inf))
  error('This functions can not take Inf as the input.');
elseif find(d ~= floor(d))
  warning('This functions is designed for converting the positive integernumbers. The results might be unexpected.');
end;

if nargin < 3
    p = 2;
end;

%***************************************************************************
% assign the length
%***************************************************************************
if nargin < 2;
    tmp = max(d);
    b1 = [];
    while tmp > 0
        b1 = [b1 rem(tmp, 2)];
        tmp = floor(tmp/2);
    end;
    n = length(b1);
end;

%*****************************************************************************
% initial value
%*****************************************************************************
b = zeros(len_d, n);

%*****************************************************************************
% parameter assignment
%*****************************************************************************
for i = 1 : len_d
    j = n;
    tmp = d(i);
    while (j >= 1) & (tmp > 0)
        b(i, j) = rem(tmp, p);
        tmp = floor(tmp/p);
        j = j - 1;
    end;
end;