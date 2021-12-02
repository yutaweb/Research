function [Datarx,Datarx_hd]=receiver(TimeSignal,ifftsize,carriers,wordsize,guardtype,...
                                     guardtime,Num_sym,Num_pilot,Doppler,proc_gain)
% TimeSignal=[1,22*80],ifftsize=64,carriers=64,wordsize=2,guardtype=2,guardtime=16,Num_sym=20,Num_pilot=2,Doppler=4,proc_gain=64
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
SymbLen=length(TimeSignal)+guardtime;% 22*80+16 ⇒ 全体の信号にもガードインターバルを付与している 	
TimeSignal=TimeSignal(1:(SymbLen-rem(SymbLen,ifftsize+guardtime)));% rem：余りを返す
% rem(22*80+16,64+16)=16 ⇒ Symblenから16を引く。 TimeSignal=22*80
numsymb=length(TimeSignal)/(ifftsize+guardtime);% 22
%***************************************************************************
% Reshape the linear time waveform into fft segments and remove guard period
%***************************************************************************
if guardtype ~= 0 % guardtypeが0でないならば、1を返し、0ならば0を返す。
   symbwaves=reshape(TimeSignal,ifftsize+guardtime,numsymb);% [80,22] 
   symbwaves=symbwaves(guardtime+1:ifftsize+guardtime,:);% [64,22] ⇒ ガードインターバルの除去
else
   symbwaves=reshape(TimeSignal,ifftsize,numsymb);
end
clear TimeSignal;%save memory

% 修正箇所
fftspect=fft(symbwaves)'./sqrt(ifftsize);% [22,64]
DataCarriers=fftspect(:,carriers);% [22,64] ⇒ FFT,IFFTは2の累乗でないとうまく適用できない
clear fftspect;%Save memory
clear symbwaves;%save memory


%*******************************
% Estimate the channel response 
%*******************************
NumCarr=size(carriers,2);% 64
H_Res=DataCarriers(1:Num_pilot,:);% [2,64] ⇒ パイロット信号:チャネル応答（伝播環境）を調べることができる h*x
if Num_pilot==1
   H_Resp=H_Res;
else
   H_Resp=sum(H_Res)/Num_pilot; % OFDMでは1ユーザーをTx:1,Rx:1としているので、伝搬環境は1つのみでOK　⇒　※1ユーザの時はチャネル推定する必要がない。
   %※MIMOの時は、伝搬環境がTx,Rxによって変化するので、注意が必要。[1,64]
   %H_Resp=DataCarriers(1,:);
end

Noise_power1=abs(DataCarriers(1,:)-H_Resp).^2;
Noise_power2=abs(DataCarriers(2,:)-H_Resp).^2;
Noise_power=(Noise_power1+Noise_power2)/2;



if Doppler==0
    data_faded1=DataCarriers(Num_pilot+1:Num_pilot+Num_sym,:);  
else
   %data_faded1=DataCarriers(Num_pilot+1:Num_pilot+Num_sym,:)./repmat(H_Resp,Num_sym,1);
   %%%ORC  
   %conj：複素共役
   data_faded1=(DataCarriers(Num_pilot+1:Num_pilot+Num_sym,:).*conj(repmat(H_Resp,Num_sym,1)))./...
               (abs(repmat(H_Resp,Num_sym,1)).^2+repmat(Noise_power,Num_sym,1)); %%MMSEC　⇒　MIMOへの拡張が難しそう。。。
   % [20,64].*[20,64]./([20,64]+[20,64])⇒要素ごとに計算を行っている。
   % 分子では、チャネル応答の自己相関を求め、チャネル応答を最大化している。分母では、分子で最大化したチャネル応答を割るとともに、ノイズの影響を軽減してあげている。
   % つまり、data_faded1には、transmitter.mでFSSシステムによって、周波数拡散したデータと軽減されたノイズが乗ったデータが入っている。
   % MMSECの方がORCよりも特性は良くなる。 ⇒ 信号を復調している。MMSEC,ORCはMLDなどと同じ役割を担っている。
   % ※MLDが最もBER特性が良いが受信アンテナの数が増加すれば、BER特性間の差は小さくなる。（マルチユーザーMIMOの教科書　参考）
end
%for xx=1:size(data_faded1,2)
%   data_faded1(:,int_pattern_1(xx))=data_faded1(:,xx);
%end



%***************************************************
% Inverse Frequenct Spreading Code(I-FSC) Operation
%***************************************************
%proc_gain=hadamard(NumCarr);
%for fsc=1:NumCarr
%   data_faded(:,fsc)=sum(data_faded1.*repmat(proc_gain(fsc,:),Num_sym,1),2)./sqrt(NumCarr);
%end

proc_gain_block=NumCarr/proc_gain;% 64/64 = 1
spread=hadamard(proc_gain);% [64,64]

if Doppler==0
    for fscb=1:proc_gain_block
      for fsc=1:proc_gain
         data_faded(:,(fscb-1)*proc_gain+fsc)=sum(data_faded1(:,(fscb-1)*proc_gain+1:fscb*proc_gain).*...
         repmat(spread(fsc,:),Num_sym,1),2)./proc_gain;
      end
   end
else
    for fscb=1:proc_gain_block% 1～1
      for fsc=1:proc_gain% 1～64
         data_faded(:,(fscb-1)*proc_gain+fsc)=sum(data_faded1(:,(fscb-1)*proc_gain+1:fscb*proc_gain).*...
         repmat(spread(fsc,:),Num_sym,1),2)./sqrt(proc_gain); % data_faded(20,1～64)=sum(data_faded1(20,64).*repmat(spread(1～64,64),20,1),2)/8⇒[20,64]
         % sum(〇,2) ⇒ 行方向の合計値を求める ⇒ 信号のみを取り出すことができる。[20,64,2]
      end
   end
end
%****************
% Data detection 
%****************
for  k = 1:Num_sym
   for cc=1:NumCarr
      if wordsize==1
         X=data_faded(k,cc);
         Xi_hd=sign(real(X));% 0より大きければ1を返し、0より小さければ-1を返す。
         Xi=real(X);% 実数部分の生データ
         Datarx_hd(k,cc)= (Xi_hd+1)/2;% 0or1のデータ
         Datarx(k,cc)= -Xi;
      elseif wordsize==2
         X=data_faded(k,cc);
         Xi_hd=sign(real(X));
         Xq_hd=sign(imag(X));
         Xi=real(X);
         Xq=imag(X);
         Datarx_hd(k,cc,1)=(Xi_hd+1)/2;
         Datarx_hd(k,cc,2)=(Xq_hd+1)/2;
         Datarx(k,cc,1)=-Xi;
         Datarx(k,cc,2)=-Xq;
      elseif wordsize==4
         X=data_faded(k,cc)*sqrt(10);
         Xi=real(X);
         Xq=imag(X);
         inst_dec1=Xi>0;
         Datarx_hd(k,cc,1)=inst_dec1;
         Datarx(k,cc,1)=-(Xi);
         if inst_dec1==0
            Xi=Xi+2; 
         else
            Xi=Xi-2; 
         end
         Datarx_hd(k,cc,2)=sign(Xi)==1;
         Datarx(k,cc,2)=-Xi;
         
         inst_dec2=Xq>0;
         Datarx_hd(k,cc,3)=inst_dec2;
         Datarx(k,cc,3)=-(Xq);
         if inst_dec2==0
            Xq=Xq+2; 
         else
            Xq=Xq-2; 
         end
         Datarx_hd(k,cc,4)=sign(Xq)==1;
         Datarx(k,cc,4)=-Xq;
      end
   end
end
clear data_faded;

