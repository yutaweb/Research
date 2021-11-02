function H_m_resp=receiver_non_FSS_MU_MIMO_AM(TimeSignal,ifftsize,carriers,wordsize,guardtype,...
                                     guardtime,Num_sym,Num_pilot,Doppler,proc_gain,Num_Tx,Num_Rx,Num_User)
% TimeSignal=[1,24*80,2,2],ifftsize=64,carriers=[1,64],wordsize=2,guardtype=2,guardtime=16,Num_pilot=2,Doppler=10,proc_gain=64,Num_User=2
% Num_Tx=4,Num_Rx=2
%****************************************************************************    
% Receiver seperates and decodes the receieved signals by using ML detection.
%
% INPUTS:
% ========
% TimeSignal=Received OFDM waveform(Tx1). 
% ifftsize=Size of ifft to use for generating the waveform
% carriers=Which carriers to use for the transmission
% wordsize=modulation command 
% guardtype=What type of guard period to use options:
%	      	0: No Guard period
%	      	1: zero level guard period
%	      	2: cyclic extension of end of symbols
% guardtime=Number of sample to use for the total guard time
% Num_sym=Number of symbols
% Num_pilot=Number of pilot signals
%
% OUTPUTS:
% ========
% Datarx1=This is the output data1 that has been decoded from the 'TimeSignal'
%
% Copyright (c) Chang-Jun Ahn 2000-2003(junny700@crl.go.jp)
%*******************************************************************************
rand('state',sum(100*clock));
%**************************************************************************
% Strip back the number of samples to make it a multiple of the symbol size
%**************************************************************************
if guardtype==0% ガードインターバルを使わない時
   guardtime=0;
end
SymbLen=zeros(1,Num_Rx,Num_User); % [0,0] ※注）zeros(2)=[[0 0];[0 0]];
for u=1:Num_User % 1～2
    for receive=1:Num_Rx % 1～2
        SymbLen(1,receive,u)=length(TimeSignal(:,:,receive,u))+guardtime;% 24*80+16 ⇒ 全体の信号にもガードインターバルを付与している
        TimeSignal(:,:,receive,u)=TimeSignal(:,1:(SymbLen(1,receive,u)-rem(SymbLen(1,receive,u),ifftsize+guardtime)),receive,u);% rem：余りを返す
        % rem(24*80+16,64+16)=16 ⇒ Symblenから16を引く。 TimeSignal=[1,24*80,2,2]
    end
end
numsymb=length(TimeSignal)/(ifftsize+guardtime);% 24



%*****************************************************************************
% Reshape the linear time waveform into fft segments and remove guard period
%*****************************************************************************
symbwaves1=zeros(ifftsize+guardtime,numsymb,Num_Rx,Num_User); % [80,24,2,2]
symbwaves=zeros(ifftsize,numsymb,Num_Rx,Num_User); % [64,24,2,4]
if guardtype ~= 0 % guardtypeが0でないならば、1を返し、0ならば0を返す。
    for u=1:Num_User
        for receive=1:Num_Rx
            symbwaves1(:,:,receive,u)=reshape(TimeSignal(:,:,receive,u),ifftsize+guardtime,numsymb);% [80,24,2,2]
            symbwaves(:,:,receive,u)=symbwaves1(guardtime+1:ifftsize+guardtime,:,receive,u);% [64,24,2,2] ⇒ ガードインターバルの除去
        end
    end
else
    for u=1:Num_User
        for receive=1:Num_Rx
            symbwaves(:,:,receive,u)=reshape(TimeSignal(:,:,receive,u),ifftsize,numsymb);
        end
    end
end
clear TimeSignal;%save memory


% 各ユーザで全データに対してFFTを掛けて、欲しいデータのみを復元し、取り出す。
fftspect=zeros(numsymb,ifftsize,Num_Rx,Num_User); % [24,64,2,2]
for u=1:Num_User
    for receive=1:Num_Rx % 1～2
        fftspect(:,:,receive,u)=fft(symbwaves(:,:,receive,u))'./sqrt(ifftsize);% [24,64,2,2]
    end
end
DataCarriers=fftspect(:,:,:,:);% [24,64,2,2] ⇒ FFT,IFFTは2の累乗でないとうまく適用できないので、fftsize=64を用いたが、計算においては冗長


clear fftspect;%Save memory
clear symbwaves;%save memory

%*******************************
% Estimate the channel response 
%*******************************
% 各ユーザからデータを送ってもらう

NumCarr=size(carriers,2);% 64 
H_Res=DataCarriers(1:Num_pilot,:,:,:);% [4,64,2,2] ⇒ パイロット信号:チャネル応答を調べることができる
% h12 → 送信アンテナ2から受信アンテナ1
% H_Res[:,:,1,1]：ユーザ1の受信アンテナ1 → h11+h12+h13+h14+noise1,h11-h12+h13-h14+noise1,h11+h12-h13-h14+noise1,h11-h12-h13+h14+noise1
% H_Res[:,:,2,1]：ユーザ1の受信アンテナ2 → h21+h22+h23+h24+noise2,h21-h22+h23-h24+noise2,h21+h22-h23-h24+noise2,h21-h22-h23+h24+noise2
% H_Res[:,:,1,2]：ユーザ2の受信アンテナ1 → h31+h32+h33+h34+noise3,h31-h32+h33-h34+noise3,h31+h32-h33-h34+noise3,h31-h32-h33+h34+noise3
% H_Res[:,:,2,2]：ユーザ2の受信アンテナ2 → h41+h42+h43+h44+noise4,h41-h42+h43-h44+noise4,h41+h42-h43-h44+noise4,h41-h42-h43+h44+noise4,
H_m1=permute(H_Res,[1 3 2 4]); % [4,2,64,2]
H_m2(:,1:2,:)=H_m1(:,:,:,1);
H_m2(:,3:4,:)=H_m1(:,:,:,2);
H_m=permute(H_m2,[2 1 3]); % [4,4,64] → 上記の順番に並んでいる（教科書の通りの並び順）
% H_m(1,:,:)：ユーザ1のアンテナ1に届いたチャネル応答
% H_m(2,:,:)：ユーザ1のアンテナ2に届いたチャネル応答
% H_m(3,:,:)：ユーザ2のアンテナ3に届いたチャネル応答
% H_m(4,:,:)：ユーザ2のアンテナ4に届いたチャネル応答
H=hadamard(Num_Tx);%[4,4]
pilot_signal_resp=repmat(H,1,1,ifftsize);

%***********************
%  channel_estimation
%***********************
% 基地局側での処理（各ユーザ間ではチャネル情報を共有できないため）→ 欲しいチャネル応答のデータが得られた
for cc=1:NumCarr
    H_m_resp(:,:,cc)=H_m(:,:,cc)*inv(pilot_signal_resp(:,:,cc));
end
% H_m_respには、以下の4×4の行列が入っている
% [h11 h12 h13 h14
%  h21 h22 h23 h24
%  h31 h32 h33 h34
%  h41 h42 h43 h44]


