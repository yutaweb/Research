function Gen_data= MIMO_data(NumSymb,wordsize,NumCarr,Num_Tx)% Numsymb=20,wordsize=1～2,NumCarr=64,Num_Tx が渡ってくる
%***********************************************************
% This function generates random data 
%
% INPUTS
% ======
% NumSymb=number of symbols
% wordsize=modulation command(1:BPSK, 2:QPSK)
% NumCarr=number of carriers
% Num_Tx=transmit antenna
% 
% OUTPUT
% ======
% Gen_data=generated bits
%
% Copyright (c) Chang-Jun Ahn 2000-2003(junny700@crl.go.jp)
%***********************************************************
rand('state',sum(100*clock));%
Gen_data=floor(rand([NumSymb,NumCarr,wordsize,Num_Tx])*2);
% randは0～1までの実数を返す。
% floorは小数点をすべて切り捨てるので、randで生成された数はこのままだと、すべて切り捨てられて0となる。
% ⇒randに2を掛けることで、[20,64,2,2]の0or1をランダムに作成している。
