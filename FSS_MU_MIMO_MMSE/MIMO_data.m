function Gen_data= MIMO_data(NumSymb,wordsize,NumCarr,Num_Tx)% Numsymb=20,wordsize=1～2,NumCarr=64,Num_Tx=4
% NumCarr=number of carriers
% Num_Tx=transmit antenna
% 
% OUTPUT
% ======
% Gen_data=generated bits
%
% Copyright (c) Chang-Jun Ahn 2000-2003(junny700@crl.go.jp)
%***********************************************************
rand('state',sum(100*clock));
% 各送信アンテナから異なるストリームを送信する場合（空間多重しているが、アンテナダイバーシチとは言わない）
Gen_data=floor(rand([NumSymb,NumCarr,wordsize,Num_Tx])*2); % [20,64,2,4]：S1, S2, S3, S4

% 各送信アンテナから同じストリームを送信する場合（アンテナダイバーシチという）
% Gen_data1=floor(rand([NumSymb,NumCarr,wordsize])*2);
% Gen_data2=floor(rand([NumSymb,NumCarr,wordsize])*2);
% Gen_data=zeros(NumSymb,NumCarr,wordsize,Num_Tx);
% 
% Gen_data(:,:,:,1)=Gen_data1;
% Gen_data(:,:,:,2)=Gen_data1;
% Gen_data(:,:,:,3)=Gen_data2;
% Gen_data(:,:,:,4)=Gen_data2;
% randは0～1までの実数を返す。
% floorは小数点をすべて切り捨てるので、randで生成された数はこのままだと、すべて切り捨てられて0となる。
% ⇒randに2を掛けることで、[20,64,2,2]の0or1をランダムに作成している。
