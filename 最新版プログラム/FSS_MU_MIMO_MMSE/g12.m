function g = g12(k)
%************************************************************
%G12 Transfer function of rate 1/2 convolutional encoder.
%   G = G12(K) outputs a binary format transfer function of rate 1/2
%   convolutional encoder.
%
%   INPUT
%   ======
%   k= constraint length
%  
%   OUTPUT
%   ========
%   g=Generator vector
%
%  Copyright (c) Chang-Jun Ahn 2000 (jun@sasase.ics.keio.ac.jp)
%************************************************************

if k == 3
    g = [1 1 1;
         1 0 1];
elseif k == 4
    g = [1 1 1 1;
         1 1 0 1];
elseif k == 5
    g = [1 1 1 0 1;
         1 0 0 1 1];
elseif k == 6
    g = [1 1 1 1 0 1;
         1 0 1 0 1 1];
elseif k == 7
    g = [1 1 1 1 0 0 1;
         1 0 1 1 0 1 1];
elseif k == 8
    g = [1 1 1 1 1 0 0 1;
         1 0 1 0 0 1 1 1];
elseif k == 9
    g = [1 1 1 1 0 1 0 1 1;
         1 0 1 1 1 0 0 0 1];
elseif k == 10
    g = [1 1 0 1 1 0 0 1 0 1;
         1 0 0 1 1 1 0 1 1 1];
elseif k == 11
    g = [1 1 1 1 0 1 1 0 0 0 1;
         1 0 0 1 1 0 1 1 1 0 1];
elseif k == 12
    g = [1 0 1 1 1 1 0 1 0 0 1 1;
         1 0 0 0 1 1 0 1 1 1 0 1];
elseif k ==13
    g = [1 1 1 1 1 1 0 1 1 0 0 0 1;
         1 0 0 0 1 0 1 0 1 1 0 1 1];
elseif k ==14
    g = [1 0 1 1 1 0 0 1 0 1 0 0 1 1;
         1 0 0 0 1 1 1 0 1 1 1 1 0 1];
else
    error('Invalid constraint length.');
end;



 