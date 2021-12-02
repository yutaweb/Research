function [OutSignal,Signal_dis, Fading]=channel_only(TxSignal,Multipath,MultiNo,Transrate,...
    Doppler,Delay,SNR, com_delay)
% TxSignal=[80,22,2],Multipath=[7,1],MultiNo=7,Transrate=1000,Doppler=4,Delay=2,SNR=7,com_delay=0,Num_Tx=2
% channel_rand=[2,2,4,7,7,rep],[Num_Tx,Num_Rx,Num_User,MultiNo,k,rep]
%*******************************************************************
% Channel applies a channel model to the waveform including 
% multipath and fading.
%
% INPUTS:
% =======
% TxSignal=Time waveform of signal to apply the channel to.  
% Multipath=Multipath magnitude if -NdB signal = 10^(0-(N/10))
% Num_Tx_ant=number of transmit antennas
% MultiNo=Number of Multipath
% Transrate=Total transmit bit rates ⇒ 1秒間に送ることのできるデータの量
% Doppler=Doppler frequency  
% Delay=Delay time
% SNR=signal to noise ratio
% wordsize=USTM command
%
% OUTPUTS:
% ==========
% TimeSignal=Output signal after the model.
%
% Copyright (c) 2000-2002 Chang-Jun Ahn (jun@sasase.ics.keio.ac.jp)
%*******************************************************************
FreqNo=size(TxSignal,1);% 80
Numsymb=size(TxSignal,2);% 22
Signal_dis=std(reshape(TxSignal,1,size(TxSignal,1)*size(TxSignal,2)));
% std：標準偏差（平均を取る時にN-1で割っていることに注意：標本なので）⇒ [1,80*22]に対して標準偏差を求めている。
%*******************
% Jakes fading model ⇒　論文で確認する
%*******************
if Doppler ~= 0% Dopplerが0でなければ、1を0であれば、0を返す。⇒ 今回はDoppler=4
   if com_delay==0
      Ns=Numsymb;% Ns=22
   else
      Ns=Numsymb*com_delay;
   end
   fs=Transrate;% fs=1000
   fm=Doppler;% fm=4
   %*************************************
   % Jakes fading model simulation start
   %*************************************
   for DD=1:MultiNo% 1～7 ⇒ マルチパス分だけフェーディングさせる
     Ts=1/Transrate;%symbol period ⇒ 1/1000   
     randpoint = rand(1);% 0～1の間をランダムにとる
     randtime = 1000*rand;
     N_0   = 32;% Number of ocsillator
     N = 4.0 * N_0;
     matrix = hadamard(N_0*2);
     tt =zeros(Ns,1);
     for t=1:Ns%
       tt(t,:)=randtime+Ts*(t-1)+tt(t,:);
       omega_M = 2.0  * pi * Doppler;
       theta_m = 2.0 * pi * randpoint;
       %***********
       % hadamard
       %***********
       for j = 1
         x_c(j) = 0.0;
         x_s(j) = 0.0;
         for  i = 1:1:N_0 % 32
           b = (pi*(i)/(N_0));
           omega_m = omega_M*cos(2.0*pi*(i-0.5)/N);
           x_c(j,i) = matrix(j,i)*(cos(b)*cos(omega_m*tt(t,:)+theta_m)...
                      +sqrt(-1)*sin(b)*cos(omega_m*tt(t,:)+theta_m));
         end
       end
        x_cc = sqrt(4.0 / N_0)*sum(x_c,2);
       fad(t,:) = [x_cc(1,1)];
       x_cc(j,1) = 0.0;
       x_ss(j,1) = 0.0;
    end
    r(DD,:)=fad';
   end%end of Jakes model
   %r=r./(repmat( (mean(r,2)),1,size(r,2)));
   if com_delay==0
      Fading=r;% [7,22]
   else
      Fading(:,1:Numsymb)=r(:,Numsymb*(com_delay-1)+1:Numsymb*com_delay);
   end
   
   
   %****************************
   % Make flat fading situation
   %****************************
   if MultiNo==1 
      fr=repmat(r(:,1:Numsymb),FreqNo,1);%Make fading replica [7,22]を80個コピー
      TxSignal=TxSignal.*fr;
      OutSignal = reshape(TxSignal,1,size(TxSignal,1)*size(TxSignal,2));% [1,80*22]
      clear TxSignal;
   else
   %******************************************************
   % Make multipath frequency selective fading situation
   %******************************************************
      Txsignal=reshape(TxSignal,1,FreqNo*Numsymb);% [1,80*22]
      clear TxSignal; 
      
      Multi=repmat(sqrt(Multipath),1,FreqNo*Numsymb);% sqrt([7,1])を[1,22*7]コピー ⇒ [7,7*22]
      fr=reshape((r(:,1:Numsymb))',1,Numsymb*MultiNo);% [1,22*7]
      fr=repmat(fr,FreqNo,1);% [80,22*7]
      fr_1=reshape(fr,1,FreqNo*Numsymb*MultiNo);% [1,80*22*7]
      clear fr;
      fr_2=(reshape(fr_1,FreqNo*Numsymb,MultiNo))';% [80*22,7]
      clear fr_1;
      Signal=repmat(Txsignal,MultiNo,1).*(Multi).*fr_2;% [80*7,22].*[7,7*22]*[80*22,7]
      clear Txsignal;
      for a=1:MultiNo
        signal(a,1:FreqNo*Numsymb)=shift(Signal(a,:),-(Delay*(a-1)));% [7,80*22]
      end
      fadcheck = signal; 
      OutSignal=sum(signal);% pathの分だけ、フェーディングさせた信号を足し合わせる ⇒ [1,80*22]
    end%end of fading situation
else% Dopplerが0の場合
    if MultiNo==1
       OutSignal=TxSignal; 
       clear TxSignal;% save the memory
       OutSignal = reshape(OutSignal,1,size(OutSignal,1)*size(OutSignal,2));
    else
       Txsignal=reshape(TxSignal,1,FreqNo*Numsymb);
       clear TxSignal;
       
       Multi=repmat(sqrt(Multipath),1,FreqNo*Numsymb);
       Signal=repmat(Txsignal,MultiNo,1).*Multi;
       clear Txsignal;
       for a=1:MultiNo
          signal(a,1:FreqNo*Numsymb)=shift(Signal(a,:),-(Delay*(a-1)));
       end
       OutSignal=sum(signal);
       clear Signal;
    end
    Fading=0;
 end
