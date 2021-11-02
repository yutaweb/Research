function [Datarx,Datarx_hd]=receiver_FSS_MIMO_OFDMA(TimeSignal,ifftsize,carriers,wordsize,guardtype,...
                                     guardtime,Num_sym,Num_pilot,Doppler,proc_gain,Num_Tx,Num_Rx,Num_User,noise,BD_H11_U1,BD_H22_U1)
% 転置：.'
% 複素共役転置：'
% TimeSignal=[1,24*64,2,2],ifftsize=64,carriers=[1,64],wordsize=1,guardtype=2,guardtime=16,Num_sym=20,
% Num_pilot=4,Doppler=10,proc_gain=64,noise=[1,64*24,2,2], H_m_resp=[4,4,64]
% Num_Tx=4,Num_Rx=2
%****************************************************************************
% Receiver seperates and decodes the receieved signals by using ML detection.
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
if guardtype==0 % ガードインターバルを使わない時
   guardtime=0;
end
SymbLen=zeros(1,Num_Rx,Num_User); % [0,0] ※注）zeros(2)=[[0 0];[0 0]];

for u=1:Num_User % 1～2
    for receive=1:Num_Rx % 1～2
        SymbLen(1,receive,u)=length(TimeSignal(:,:,receive,u))+guardtime;% 22*80+16 ⇒ 全体の信号にもガードインターバルを付与している
        TimeSignal(:,:,receive,u)=TimeSignal(:,1:(SymbLen(1,receive,u)-rem(SymbLen(1,receive,u),ifftsize+guardtime)),receive,u);% rem：余りを返す
        % rem(22*80+16,64+16)=16 ⇒ Symblenから16を引く。 TimeSignal=[1,24*80,2,4]
    end
end

numsymb=length(TimeSignal)/(ifftsize+guardtime);% 24

%*****************************************************************************
% Reshape the linear time waveform into fft segments and remove guard period
%*****************************************************************************
symbwaves1=zeros(ifftsize+guardtime,numsymb,Num_Rx,Num_User); % [80,24,2,2]
symbwaves=zeros(ifftsize,numsymb,Num_Rx,Num_User); % [64,24,2,2]
if guardtype ~= 0 % guardtypeが0でないならば、1を返し、0ならば0を返す。
    for u=1:Num_User
        for receive=1:Num_Rx
            symbwaves1(:,:,receive,u)=reshape(TimeSignal(:,:,receive,u),ifftsize+guardtime,numsymb);% [80,24,2,2]
            symbwaves(:,:,receive,u)=symbwaves1(guardtime+1:ifftsize+guardtime,:,receive,u);% [64,24,2,2]
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


fftspect=zeros(numsymb,ifftsize,Num_Rx,Num_User); % [24,64,2,2] ※パイロット信号部分のみが一致する可能性がある
for u=1:Num_User
    for receive=1:Num_Rx
        fftspect(:,:,receive,u)=fft(symbwaves(:,:,receive,u))'./sqrt(ifftsize);% [24,64,2,2]
    end    
end
DataCarriers=fftspect(:,:,:,:);% [24,64,2,2] ⇒ FFT,IFFTは2の累乗でないとうまく適用できないので、fftsize=64を用いたが、計算においては冗長

clear fftspect;%Save memory
clear symbwaves;%Save memory

%***********************************
%　　　　　 beamforming
%***********************************
% BD-SVD法
% DataCarriers1=zeros(numsymb,ifftsize,Num_Rx,Num_User); % [24,64,2,2]
% for symbol=1:numsymb
%     for cc=1:ifftsize
%         DataCarriers1(symbol,cc,:,1)=BD_H11_U1(:,:,cc)'*permute(DataCarriers(symbol,cc,:,1),[3,1,2,4]);
%         DataCarriers1(symbol,cc,:,2)=BD_H22_U1(:,:,cc)'*permute(DataCarriers(symbol,cc,:,2),[3,1,2,4]);
%     end
% end

% CI法（特に操作は不要）
DataCarriers1=DataCarriers;
%***********************************
%　　　　　 beamforming
%***********************************

%[2,2,64,4]⇒各ユーザで独立したチャネル推定値が入っている。

proc_gain_block=ifftsize/proc_gain;% 64/64 = 1
data_faded1=DataCarriers1(Num_pilot+1:Num_pilot+Num_sym,:,:,:); % [20,64,2,2]
data_faded=zeros(Num_sym,ifftsize,Num_Rx,Num_User); % [20,64,2,2]
data_faded2=zeros(Num_sym,ifftsize,Num_Tx); % [20,64,4]


%***************************************************
% Inverse Frequenct Spreading Code(I-FSC) Operation
%***************************************************
spread=hadamard(proc_gain);% [64,64]

if Doppler==0
    for u=1:Num_User
        for i=1:Num_Rx
            for fscb=1:proc_gain_block
                for fsc=1:proc_gain
                    data_faded(:,(fscb-1)*proc_gain+fsc,i,u)=sum(data_faded1(:,(fscb-1)*proc_gain+1:fscb*proc_gain.i,u).*...
                        repmat(spread(fsc,:),Num_sym,1),2)./proc_gain;
                end
            end
        end
    end
else
    for u=1:Num_User
        for i=1:Num_Rx
            for fscb=1:proc_gain_block% 1～1
                for fsc=1:proc_gain% 1～64
                    data_faded(:,(fscb-1)*proc_gain+fsc,i,u)=sum(data_faded1(:,(fscb-1)*proc_gain+1:fscb*proc_gain,i,u).*...
                        repmat(spread(fsc,:),Num_sym,1),2)./sqrt(proc_gain); % data_faded(20,1～64,1～2,1～4)=sum(data_faded1(20,1～64,1～2,1～4).*repmat(spread(1～64,64),20,1),2)/8⇒[20,64]
                    % sum(〇,2) ⇒ 行方向の合計値を求める ⇒ 信号のみを取り出すことができる。[20,64,2,4]
                    % 夫々のユーザは独立しているので、各ユーザでFSS逆拡散をする
                    % [20,64,2,2]
                end
            end
        end
    end
end
data_faded2(:,:,1:2)=data_faded(:,:,:,1);
data_faded2(:,:,3:4)=data_faded(:,:,:,2);

%****************
% Data detection 
%****************
for i=1:Num_Tx % 1 ～ 4
% for i=1:1 % 1 ～ 4
    for k = 1:Num_sym % 1 ～ 20
       for cc=1:ifftsize % 1 ～ 64
          if wordsize==1
             X=data_faded2(k,cc);
             Xi_hd=sign(real(X));% 0より大きければ1を返し、0より小さければ-1を返す。
             Xi=real(X);% 実数部分の生データ
             Datarx_hd(k,cc)= (Xi_hd+1)/2;% 0or1のデータ
             Datarx(k,cc)= -Xi;
          elseif wordsize==2
             % 送信アンテナでは、Dataは[NumSymb,NumCarr,wordsize,Num_Tx]となっている。（-1or1のデータ）
             % 以下の復調データの格納がうまくいっていない。
             X(k,cc,i)=data_faded2(k,cc,i)*sqrt(2); % [20,64,1～4]
             Xi_hd(k,cc,i)=sign(real(X(k,cc,i))); % [20,64,1～4]
             Xq_hd(k,cc,i)=sign(imag(X(k,cc,i))); % [20,64,1～4]
             Xi(k,cc,i)=real(X(k,cc,i)); % [20,64,1～4]
             Xq(k,cc,i)=imag(X(k,cc,i)); % [20,64,1～4]
             Datarx_hd(k,cc,1,i)=(Xi_hd(k,cc,i)+1)/2; 
             Datarx_hd(k,cc,2,i)=(Xq_hd(k,cc,i)+1)/2;
             Datarx(k,cc,1,i)=-Xi(k,cc,i);
             Datarx(k,cc,2,i)=-Xq(k,cc,i);
          elseif wordsize==4
             X=data_faded2(k,cc)*sqrt(5);
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

