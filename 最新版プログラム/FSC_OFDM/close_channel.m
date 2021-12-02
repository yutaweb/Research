function [OutSignal,Signal_dis]=close_channel(TxSignal,Multipath,MultiNo,Delay,r,Doppler);
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
% Transrate=Total transmit bit rates
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
FreqNo=size(TxSignal,1);
Numsymb=size(TxSignal,2);
Signal_dis=std(reshape(TxSignal,1,size(TxSignal,1)*size(TxSignal,2)));
rand('state',sum(100*clock));
%*******************
%Jakes fading model
%*******************
   
if Doppler ~= 0;
   %****************************
   % Make flat fading situation
   %****************************
   if MultiNo==1 
      fr=repmat(r,FreqNo,1);%Make fading replica 
      TxSignal=TxSignal.*fr;
      OutSignal = reshape(TxSignal,1,size(TxSignal,1)*size(TxSignal,2));
      clear TxSignal;
   else
   %******************************************************
   % Make multipath frequency selective fading situation
   %******************************************************
      Txsignal=reshape(TxSignal,1,FreqNo*Numsymb);
      clear TxSignal; 
      
      Multi=repmat(sqrt(Multipath),1,FreqNo*Numsymb);
      fr=reshape((r)',1,Numsymb*MultiNo);
      fr=repmat(fr,FreqNo,1);
      fr_1=reshape(fr,1,FreqNo*Numsymb*MultiNo);
      clear fr;
      fr_2=(reshape(fr_1,FreqNo*Numsymb,MultiNo))';
      clear fr_1;
      Signal=repmat(Txsignal,MultiNo,1).*(Multi).*fr_2;
      clear Txsignal;
      for a=1:MultiNo
        signal(a,1:FreqNo*Numsymb)=shift(Signal(a,:),-(Delay*(a-1)));
      end
      fadcheck = signal; 
      OutSignal=sum(signal);
    end%end of fading situation
 else
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
 end
