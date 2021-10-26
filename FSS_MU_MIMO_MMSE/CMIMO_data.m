function [Gen_data,Ori_data,int_pattern]=MIMO_data(NumSymb,wordsize,NumCarr,g,Num_Tx,Num_User)
% NumSymb=20,wordsize=2（QPSK）,NumCarr=64,g=0or1 [2,7],Num_Tx=4,Num_user=2
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

%% 全ストリームに対して一括でFEC及びインタリーブを掛ける場合（各ユーザごとにかける）
Total_bit=(NumSymb*NumCarr*wordsize*Num_Tx/Num_User)/2-mem;% 20*64*2-6
Ori_data=randi([0,1],1,Total_bit); % [1,20*64*2-6]
Data=zeros(1,NumSymb*NumCarr*wordsize*Num_Tx/Num_User,Num_User); % [1,20*64*2*2,2]
temp_1=zeros(NumSymb*NumCarr*wordsize*Num_Tx/Num_User,Num_User); % [20*64*2*2,2]
int_pattern=zeros(NumSymb*NumCarr*wordsize*Num_Tx/Num_User,Num_User); %[20*64*2*2,2]
for u=1:Num_User
    Data(:,:,u) = encoder(g,Ori_data); % g=[2,7],Ori_data=[1,20*64*2-6] 返却データ量は、[1,20*64*2*2]
    [temp_1(:,u), int_pattern(:,u)]=sort(rand(1,size(Data,2)));% temp_1：値,int_pattern：インデックス
    for xx=1:size(Data,2) % 1～20*64*2*2
       Datatx(1,xx,u)=Data(1,int_pattern(xx,u),u); % ランダムに並べ替えている
    end
    Datatx1(:,:,u)=reshape(Datatx(1,:,u),[NumSymb,NumCarr*wordsize*Num_Tx/Num_User]);% [20,64*2*2]
    if wordsize==1
       Gen_data=Datatx1;
    else
       Gen_data1(:,:,:,u)=reshape(Datatx1(:,:,u),[NumSymb,NumCarr,wordsize*Num_Tx/Num_User]);% [20,64,2*2]
       Gen_data2(:,:,:,:,u)=reshape(Gen_data1(:,:,:,u),[NumSymb,NumCarr,wordsize,Num_Tx/Num_User]); % [20,64,2,2]
    end
end

if wordsize ~= 1
     Gen_data(:,:,:,1:2)=Gen_data2(:,:,:,:,1);
     Gen_data(:,:,:,3:4)=Gen_data2(:,:,:,:,2);
end


% Gen_dataについての説明
% → [データシンボル数,サブキャリア数,変調方式,ユーザ毎の送信アンテナ数,ユーザ数]

clear Datatx;%save memory
clear data;%save memory
