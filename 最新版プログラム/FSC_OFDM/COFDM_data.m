function [Gen_data,Ori_data,int_pattern]=OFDM_data(NumSymb,wordsize,NumCarr,g)
% NumSymb=20,wordsize=1～2（QPSK）,NUmCarr=64,g=0or1
%***********************************************************
% This function generates random data 
%
% INPUTS
% ======
% NumSymb=number of symbols
% wordsize=modulation command(1:BPSK, 2:QPSK)
% NumCarr=number of carriers
% g=generate matrix
%
% OUTPUT
% ======
% Gen_data=generated bits
% Ori_data=Generated original bits
%
% Copyright (c) Chang-Jun Ahn 2000-2003(junny700@crl.go.jp)
%***********************************************************
rand('state',sum(100*clock));%
mem=size(g,2)-1;% strength=7の時、[2,7] size(g,2) ⇒ 7 よって、6
Total_bit=((NumSymb*NumCarr*wordsize)/2)-mem;% 20*64-6
% Ori_data=randint(1,Total_bit); 【修正済み】
Ori_data=randi([0,1],1,Total_bit);% [1,20*64-6]の0or1
Data = encoder(g,Ori_data);% g=[2,7],Ori_data=[1,20*64-6]　⇒ 20*64*2（2倍のデータ量が返ってくる）
[temp_1, int_pattern]=sort(rand(1,size(Data,2)));% temp_1:値,int_pattern:index ⇒　データをランダムに並べ替えることで、雑音に強くなる。
for xx=1:size(Data,2)% 1～20*64*2
   Datatx(1,xx)=Data(1,int_pattern(xx)); % データがランダムに並べ替えられている。
end
Datatx=reshape(Datatx,NumSymb,NumCarr*wordsize); % [20,64*2] 
if wordsize==1
   Gen_data=Datatx;
else
   Gen_data=reshape(Datatx,[NumSymb,NumCarr,wordsize]); % [20,64,2]
end
clear Datatx;%save memory
clear data;%save memory
