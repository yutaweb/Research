function word=receiver(TimeSignal,ifftsize,carriers,wordsize,guardtype,...
                         guardtime,Num_sym,Num_pilot,Doppler,Num_Rx,Num_Tx,noise)
% TimeSignal=[1,80*22],ifftsize=64,carriers=64,wordsize=1～2,guardtype=2,guardtime=16,Num_sym=20,Num_pilot=2,Doppler=10
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
%Strip back the number of samples to make it a multiple of the symbol size
%**************************************************************************
if guardtype==0
   guardtime=0;
end
SymbLen=zeros(1,Num_Rx); % [0,0] ※注）zeros(2)=[[0 0];[0 0]];
for receive=1:Num_Rx % 1～2
    SymbLen(1,receive)=length(TimeSignal(:,:,receive))+guardtime;% 22*80+16 ⇒ 全体の信号にもガードインターバルを付与している 	
    TimeSignal(:,:,receive)=TimeSignal(:,1:(SymbLen(1,receive)-rem(SymbLen(1,receive),ifftsize+guardtime)),receive);% rem：余りを返す
    % rem(22*80+16,64+16)=16 ⇒ Symblenから16を引く。 TimeSignal=[1,22*80,2]
end
numsymb=length(TimeSignal)/(ifftsize+guardtime);% 22

%***************************************************************************
%Reshape the linear time waveform into fft segments and remove guard period
%***************************************************************************
symbwaves1=zeros(ifftsize+guardtime,numsymb,Num_Rx); % [80,22,2]
symbwaves=zeros(ifftsize,numsymb,Num_Rx); % [64,22,2]
if guardtype ~= 0 % guardtypeが0でないならば、1を返し、0ならば0を返す。
    for receive=1:Num_Rx
        symbwaves1(:,:,receive)=reshape(TimeSignal(:,:,receive),ifftsize+guardtime,numsymb);% [80,22,2] 
        symbwaves(:,:,receive)=symbwaves1(guardtime+1:ifftsize+guardtime,:,receive);% [64,22,2] ⇒ ガードインターバルの除去
    end
else
    for receive=1:Num_Rx
        symbwaves(:,:,receive)=reshape(TimeSignal(:,:,receive),ifftsize,numsymb);
    end
end
clear TimeSignal;%save memory

fftspect=zeros(numsymb,ifftsize,Num_Rx); % [22,64,2]
for receive=1:Num_Rx
    fftspect(:,:,receive)=fft(symbwaves(:,:,receive))'./sqrt(ifftsize);% [22,64,2]
end
DataCarriers=fftspect(:,:,:);% [22,64]
clear fftspect;%Save memory
clear symbwaves;%save memory


%*******************************
% Estimate the channel response 
%*******************************
NumCarr=size(carriers,2);% 64 
H_Res=DataCarriers(1:Num_pilot,:,:);% [2,64,2] ⇒ パイロット信号:チャネル応答を調べることができる
h1_1=permute(H_Res(1,:,1),[1 3 2]);% h1+h2'+n1
h1_2=permute(H_Res(2,:,1),[1 3 2]);% h1-h2'+n1'
h2_1=permute(H_Res(1,:,2),[1 3 2]);% h1'+h2+n2
h2_2=permute(H_Res(2,:,2),[1 3 2]);% h1'-h2+n2 '[1,1,64]
H_m=[h1_1 h1_2;h2_1 h2_2];% [2,2,64] ← MIMOにおいてはパス間の相関が全然ないので、


H=hadamard(Num_Tx);%[2,2]
pilot_signal_resp=repmat(H,1,1,ifftsize);
for cc=1:NumCarr
    H_m_resp(:,:,cc)=H_m(:,:,cc)*inv(pilot_signal_resp(:,:,cc));
end

gause_noise=Num_Tx*(std(noise(:,:,1)))^2.*eye(Num_Rx);
Estimated_SNR=10*log10(gause_noise); % [1,,64]
disp(Estimated_SNR);
if wordsize==4
    Average_SNR=Estimated_SNR-10*log10(4)-10*log10(2);
else
    Average_SNR=Estimated_SNR-10*log10(4); % 平均をSNRを求めている
end


if Average_SNR>=11
   word=4;
else
   word=2;
end

