function y = encoder(g, x)
% g=[7,2],x=Data [1,20*64*2-6]
%****************************************************************
% Usage: y = encoder(g, x)
% This function takes as input an entire block of information bits 'x'
% (which are arranged in a row vector), and the coeficients of the 
% generator polynomials 'g', and returns as output an entire 
% convolutionally encoded codeword 'y'.  Tail bits are automatically 
% appended to force the encoder back to the all-zeros state.
%
% The encoder matrix 'g' has n rows and K columns, where the code rate
% is 1/n and the constraint length is K.  Each row corresponds to one
% of the n encoder polynomials and each column of the row is a 
% coeffiecient of the encoder polynomial.
%
% Copyright (c) Chang-Jun Ahn 2000 (jun@sasase.ics.keio.ac.jp)
%******************************************************************

[n,K] = size(g); % n=2,K=7
m = K - 1; % 6
[temp,L_info] = size(x); % temp=1,L_info=20*64*2-6

%******************************************************************
% initialize the state vector
%******************************************************************
state = zeros(1,m); % [1,6]

%******************************************************************
% zero pad the codeword
%******************************************************************
x = [x zeros(1,m)]; % [1,20*64*2] ⇒ 後ろの6個のデータを0で埋めた
L_total = L_info+m; % 20*64*2

%******************************************************************
% generate the codeword
%******************************************************************
for i = 1:L_total % 1 ～ 20*64*2
   input_bit = x(1,i);
   [output_bits, state] = encode_bit(g, input_bit, state);
   y(n*(i-1)+1:n*i) = output_bits;
end
% cording_rate=1/2なので出力値は2倍になって出てくる。


