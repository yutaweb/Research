function [Gen_data,Ori_data,int_pattern]=MIMO_data(NumSymb,wordsize,NumCarr,g,Num_Tx)
% NumSymb=20,wordsize=2（QPSK）,NumCarr=64,g=0or1 [2,7],Num_Tx=2
%***********************************************************
% This function generates random data 
%
% INPUTS
% ======
% NumSymb=number of symbols
% wordsize=modulation command(1:BPSK, 2:QPSK)
% NumCarr=number of carriers
% g=generate matrix
% Num_Tx=transmit antenna
%
% OUTPUT
% ======
% Gen_data=generated bits
% Ori_data=Generated original bits
%
% Copyright (c) Chang-Jun Ahn 2000-2003(junny700@crl.go.jp)
%***********************************************************
rand('state',sum(100*clock));
mem=size(g,2)-1;% strength=7の時、[2,7] size(g,2) ⇒ 7 よって、6
Total_bit=(NumSymb*NumCarr*wordsize*Num_Tx)/2-mem;% 20*64*2-6
% Ori_data=randint(1,Total_bit); 【修正済み】
Ori_data=randi([0,1],1,Total_bit); % [1,20*64*2-6]
Data = encoder(g,Ori_data); % g=[2,7],Ori_data=[1,20*64*2-6] 返却データ量は、[1,20*64*2*2]
[temp_1, int_pattern]=sort(rand(1,size(Data,2)));% temp_1：値,int_pattern：インデックス
for xx=1:size(Data,2) % 1～20*64*2*2
   Datatx(1,xx)=Data(1,int_pattern(xx)); % ランダムに並べ替えている
end
Datatx1=reshape(Datatx,[NumSymb,NumCarr*wordsize,Num_Tx]);% [20,64*2,2]
if wordsize==1
   Gen_data=Datatx1;
else
   Gen_data=reshape(Datatx1,[NumSymb,NumCarr,wordsize,Num_Tx]);% [20,64,2,2]
end
clear Datatx;%save memory
clear data;%save memory
