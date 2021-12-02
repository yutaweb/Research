function [Datarx,Datarx_hd]=receiver(TimeSignal,ifftsize,carriers,wordsize,guardtype,...
                                     guardtime,Num_sym,Num_pilot,Doppler,proc_gain,Num_Tx,Num_Rx,noise)
% TimeSignal=[1,22*64,2],ifftsize=64,carriers=[1,64],wordsize=2,guardtype=2,guardtime=4,Num_sym=20,Num_pilot=2,Doppler=4,proc_gain=64
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
SymbLen=zeros(1,Num_Rx); % [0,0] ※注）zeros(2)=[[0 0];[0 0]];
for receive=1:Num_Rx % 1～2
    SymbLen(1,receive)=length(TimeSignal(:,:,receive))+guardtime;% 22*80+16 ⇒ 全体の信号にもガードインターバルを付与している 	
    TimeSignal(:,:,receive)=TimeSignal(:,1:(SymbLen(1,receive)-rem(SymbLen(1,receive),ifftsize+guardtime)),receive);% rem：余りを返す
    % rem(22*80+16,64+16)=16 ⇒ Symblenから16を引く。 TimeSignal=[1,22*80,2]
end
numsymb=length(TimeSignal)/(ifftsize+guardtime);% 22

%*****************************************************************************
% Reshape the linear time waveform into fft segments and remove guard period
%*****************************************************************************
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
DataCarriers=fftspect(:,:,:);% [22,64,2] ⇒ FFT,IFFTは2の累乗でないとうまく適用できないので、fftsize=64を用いたが、計算においては冗長

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

%[2,2,64]⇒チャネル推定値が入っている。
if Num_pilot==1
   H_Resp=H_Res;
else
   H_Resp=sum(H_m_resp,[2]); % [2,1,64]
   %OFDMの場合は伝搬チャネルが1つだけだったので容易に考えることができた。
   %H_Resp=DataCarriers(1,:);
end

gause_noise=Num_Tx*(std(noise))^2.*eye(Num_Rx);
if Doppler==0
    data_faded1=DataCarriers(Num_pilot+1:Num_pilot+Num_sym,:); % [20,64,2]  
else
   
   DataCarriers1=permute(DataCarriers(Num_pilot+1:Num_pilot+Num_sym,:,:),[3 1 2]);%[2 20 64];
   % ココでの計算は、data_faded1が[20,64,2]の状態を保ったまま 
   data_faded1=zeros(Num_sym,NumCarr,Num_Rx);
   for k=1:Num_sym % 1 ～ 20
        for cc=1:NumCarr % 1 ～ 64
            data_faded1(k,cc,:)=inv(H_m_resp(:,:,cc)'*H_m_resp(:,:,cc)+gause_noise)*H_m_resp(:,:,cc)'*DataCarriers1(:,k,cc);
        end
    end
end

%****************
% Data detection 
%****************
for i=1:Num_Rx
    for  k = 1:Num_sym % 1:20
       for cc=1:NumCarr % 1:64
          if wordsize==1
             X=data_faded1(k,cc);
             Xi_hd=sign(real(X));% 0より大きければ1を返し、0より小さければ-1を返す。
             Xi=real(X);% 実数部分の生データ
             Datarx_hd(k,cc)= (Xi_hd+1)/2;% 0or1のデータ
             Datarx(k,cc)= -Xi;
          elseif wordsize==2
             % 送信アンテナでは、Dataは[NumSymb,NumCarr,wordsize,Num_Tx]となっている。（-1or1のデータ）
             % 以下の復調データの格納がうまくいっていない。
             X(k,cc,i)=data_faded1(k,cc,i); % [20,64,1～2]
             Xi_hd(k,cc,i)=sign(real(X(k,cc,i))); % [20,64,1～2]
             Xq_hd(k,cc,i)=sign(imag(X(k,cc,i))); % [20,64,1～2]
             Xi(k,cc,i)=real(X(k,cc,i)); % [20,64,1～2]
             Xq(k,cc,i)=imag(X(k,cc,i)); % [20,64,1～2]
             Datarx_hd(k,cc,1,i)=(Xi_hd(k,cc,i)+1)/2; 
             Datarx_hd(k,cc,2,i)=(Xq_hd(k,cc,i)+1)/2;
             Datarx(k,cc,1,i)=-Xi(k,cc,i);
             Datarx(k,cc,2,i)=-Xq(k,cc,i);
          elseif wordsize==4
             X=data_faded1(k,cc)*sqrt(5);
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
end
clear data_faded;

