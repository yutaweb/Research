function bin_state = bin_state( int_state, m )
% Usage: bin_state = bin_state( int_state, m )
%
% Converts an integer 'int_state' into a row vector 
% 'bin_state' of m bits.
%
% Original author: M.C. Valenti
% For academic use only

for i = m:-1:1
   state(m-i+1) = fix( int_state/ (2^(i-1)) );
   int_state = int_state - state(m-i+1)*2^(i-1);
end

bin_state = state;

