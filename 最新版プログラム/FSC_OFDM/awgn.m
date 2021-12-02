function TimeSignal=awgn(OutSignal,Signal_dis,MultiNo,Doppler,SNR,Num_sym,...
                         NumCarr,guardtime,Num_pilot,coding_rate,perfect,wordsize)
%*******************************************************************************         
% This function generates a Gaussian noise. 
%
% INPUTS:
% ========
% OutSignal=Input data.
% Signal_dis=Standard deviation of the received signals.
% MultiNo=Number of multipaths
% Doppler=Doppler frequency
% SNR=Signal to noise ratio
% Num_sym=Number of symbols
% NumCarr=Number of carriers
% guardtime=Guard time
% Num_pilot=Number of pilots
% coding_rate=coding rate
% perfect=perfect channel estimation option
% wordsize=modulation command
%			
% OUTPUTS:
% ========
% TimeSignal=Noised OFDM signal. 
%	  
% Copyright (c) 2000-2003 Chang-Jun Ahn (junny700@crl.go.jp)
%*******************************************************************************         

%**********
% add AWGN 
%**********
randn('state',sum(100*clock));
rand('state',sum(100*clock));
Time_length=length(OutSignal);%length of time signal ⇒ [1,22*64]

if wordsize~=4
    SNR=SNR+10*log10(coding_rate);% 今回はcoding_rate=1
elseif wordsize==4
    SNR=SNR+10*log10(coding_rate)+10*log10(2);%
end
if Doppler==0
   SigPow = std(OutSignal);
   NoiseFac=10^(-SNR/20);%
else
   if MultiNo==1
      SigPow = Signal_dis;
      NoiseFac = 10^(-SNR/20);
   else
      SigPow = std(OutSignal);
      NoiseFac = 10^(-SNR/20);
   end
end
if perfect==1
   noise=zeros(1,length(OutSignal));
   %noise(length(OutSignal)*(Num_pilot/(Num_sym+Num_pilot))+1:length(OutSignal))=sqrt(1/2)*(randn(1,length(OutSignal)*(Num_sym/(Num_sym+Num_pilot)))*(i)+randn(1,length(OutSignal)*(Num_sym/(Num_sym+Num_pilot)))).*(NoiseFac);%
   noise(length(OutSignal)*((Num_pilot-1)/(Num_sym+Num_pilot))+1:length(OutSignal))=sqrt(1/2)*(randn(1,length(OutSignal)*((Num_sym+1)/(Num_sym+Num_pilot)))*(i)+randn(1,length(OutSignal)*((Num_sym+1)/(Num_sym+Num_pilot)))).*(NoiseFac);%
else
   noise=sqrt(1/2)*(randn(1,length(OutSignal))*(i)+randn(1,length(OutSignal))).*(NoiseFac);%
end
TimeSignal=OutSignal(1:Time_length)+noise;
clear OutSignal;%clear SigPow;%save memory
TimeSignal = reshape(TimeSignal,1,size(TimeSignal,1)*size(TimeSignal,2));

