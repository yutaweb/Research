function H_m_resp=receiver(TimeSignal,ifftsize,carriers,wordsize,guardtype,...
                                     guardtime,Num_sym,Num_pilot,Doppler,proc_gain,Num_Tx,Num_Rx,Num_User)
% TimeSignal=[1,22*80,2,4],ifftsize=64,carriers=[1,64],wordsize=2,guardtype=2,guardtime=16,Num_pilot=2,Doppler=10,proc_gain=64,Num_User=4
% Num_Tx=2,Num_Rx=2
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
for u=1:Num_User % 1～4
    for receive=1:Num_Rx % 1～2
        SymbLen(1,receive,u)=length(TimeSignal(:,:,receive,u))+guardtime;% 28*80+16 ⇒ 全体の信号にもガードインターバルを付与している
        TimeSignal(:,:,receive,u)=TimeSignal(:,1:(SymbLen(1,receive,u)-rem(SymbLen(1,receive,u),ifftsize+guardtime)),receive,u);% rem：余りを返す
        % rem(8*80+16,64+16)=16 ⇒ Symblenから16を引く。 TimeSignal=[1,28*80,2,4]
    end
end
numsymb=length(TimeSignal)/(ifftsize+guardtime);% 22

%*****************************************************************************
% Reshape the linear time waveform into fft segments and remove guard period
%*****************************************************************************
symbwaves1=zeros(ifftsize+guardtime,numsymb,Num_Rx,Num_User); % [80,22,2,4]
symbwaves=zeros(ifftsize,numsymb,Num_Rx); % [64,22,2,4]
if guardtype ~= 0 % guardtypeが0でないならば、1を返し、0ならば0を返す。
    for u=1:Num_User
        for receive=1:Num_Rx
            symbwaves1(:,:,receive,u)=reshape(TimeSignal(:,:,receive,u),ifftsize+guardtime,numsymb);% [80,22,2,4]
            symbwaves(:,:,receive,u)=symbwaves1(guardtime+1:ifftsize+guardtime,:,receive,u);% [64,22,2,4] ⇒ ガードインターバルの除去
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

% 各ユーザで全データに対してFFTを掛けて、欲しいでたのみを復元し、取り出す。
fftspect=zeros(numsymb,ifftsize,Num_Rx,Num_User); % [22,64,2,4]
for u=1:Num_User
    for receive=1:Num_Rx % 1～2
        fftspect(:,:,receive,u)=fft(symbwaves(:,:,receive,u))'./sqrt(ifftsize);% [22,64,2,4]
    end
end
DataCarriers=fftspect(:,:,:,:);% [22,64,2,4] ⇒ FFT,IFFTは2の累乗でないとうまく適用できないので、fftsize=64を用いたが、計算においては冗長

clear fftspect;%Save memory
clear symbwaves;%save memory

%*******************************
% Estimate the channel response 
%*******************************
% 各ユーザからh1_1,h1_2といったデータを送ってもらう
NumCarr=size(carriers,2);% 64 
H_Res=DataCarriers(1:Num_pilot,:,:,:);% [2,64,2,4] ⇒ パイロット信号:チャネル応答を調べることができる
H_m=permute(H_Res,[3 1 2 4]); % [2,2,64,4]
H=hadamard(Num_Tx);%[2,2]
pilot_signal_resp=repmat(H,1,1,ifftsize);
for u=1:Num_User
    for cc=1:NumCarr
        H_m_resp(:,:,cc,u)=H_m(:,:,cc,u)*inv(pilot_signal_resp(:,:,cc));
    end
end



