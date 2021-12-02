function word=receiver(TimeSignal,ifftsize,carriers,wordsize,guardtype,...
                         guardtime,Num_sym,Num_pilot,Doppler)
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
SymbLen=length(TimeSignal)+guardtime; % 22*80+16
TimeSignal=TimeSignal(1:(SymbLen-rem(SymbLen,ifftsize+guardtime))); % 22*80
numsymb=length(TimeSignal)/(ifftsize+guardtime); % 22

%***************************************************************************
%Reshape the linear time waveform into fft segments and remove guard period
%***************************************************************************
if guardtype ~= 0 % guardtypeが0でなければ、1でそうでなければ、0
   symbwaves=reshape(TimeSignal,ifftsize+guardtime,numsymb); % [80,22]
   symbwaves=symbwaves(guardtime+1:ifftsize+guardtime,:);% [64,22]
else
   symbwaves=reshape(TimeSignal,ifftsize,numsymb);
end
clear TimeSignal;%save memory

fftspect=fft(symbwaves)';% [22,64]
DataCarriers=fftspect(:,carriers);% [22,64]
clear fftspect;%Save memory
clear symbwaves;%save memory


%*******************************
% Estimate the channel response 
%*******************************
NumCarr=size(carriers,2); % 64
H_Res=DataCarriers(1:Num_pilot,:); % [2,64]
if Num_pilot==1
   H_Resp=H_Res;
else
   H_Resp=sum(H_Res)/Num_pilot; % [1,64] ⇒ パイロット信号の平均を取った（チャネル応答）
end

Noise_power1=abs(DataCarriers(1,:)-H_Resp).^2;
Noise_power2=abs(DataCarriers(2,:)-H_Resp).^2;
Noise_power=(Noise_power1+Noise_power2)/2;
Estimated_SNR=10*log10(abs(H_Resp).^2./Noise_power); % [1,,64]
if wordsize==4
    Average_SNR=sum(Estimated_SNR,2)/size(Estimated_SNR,2)-10*log10(4)-10*log10(2);
else
    Average_SNR=sum(Estimated_SNR,2)/size(Estimated_SNR,2)-10*log10(4); % 平均をSNRを求めている
end


if Average_SNR>=11
   word=4;
else
   word=2;
end

