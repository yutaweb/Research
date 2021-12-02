function [output, state] = encode_bit(g, input, state)
% g=[2,7].input=1,state=[1,6]
% Usage: [output, state] = encode_bit(g, input, state)
%
% Takes as an input a single bit to be encoded 'input', 
% the coeficients of the generator polynomials 'g', and
% the current state vector 'state'.
% Returns as output n encoded data bits 'output' (where 1/n 
% is the code rate), and the new state vector 'state'.
%

% the rate is 1/n
% k is the constraint length
% m is the amount of memory
[n,k] = size(g); % n=2,k=7
m = k-1; % 6

% determine the next output bit
for i=1:n % n=2
   output(i) = g(i,1)*input;
   for j = 2:k
      output(i) = xor(output(i),g(i,j)*state(j-1)); 
      %xorによって、0or1を生成している。 ⇒ 互いに異なる時のみ1を取る。
      %g=[1 1 1 1 0 0 1;
      %   1 0 1 1 0 1 1]
   end
end

state = [input, state(1:m-1)]; % [input_bit,0,0,0,0,0]　みたいな感じで更新されていく。

