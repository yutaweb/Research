%******************************************************************
% This program simulates an OFDM with multi-path signals.
% Multi-paths consist of multiple reflection.
% This program considers BPSK,QPSK, and 16QAM modulations.
% This OFDM can obtain the maximum throughput of 10Mbps(16QAM)
%
% Modification
% =============
% First Modification: 23/06/2000(OFDM)
% Second Modification:24/07/2000(BPSK, QPSK)
% Third Modification:14/08/2000(16QAM)
% Fourth Modification:21/08/2002
% Fifth Modification:11/02/2005(FSC)
% 
% Copyright (c) 2000-2005 Chang-Jun Ahn (junny700@nict.go.jp)
%******************************************************************

%******************
% Simulation start
%******************
clear all;
% flops(0);
tic;%Measure the time it takes to run the simulation
rand('state',sum(100*clock));%random probability

%***************************
% Definition for simulation
%***************************
Doppler=10;
com_delay=0;

ifftsize=64;
guardtype=2;
guardtime=16;
windowtype=0;
CarrSpacing=1; % ここを変更することで、carriersのサイズを変更することができる。

MultiNo=7;
Transrate=78100*2;
Delay=2;%Number of delay sample
ww=2;%dB degradation

NumCarr=64;   
Num_sym=20;%Number of data symbols    
Num_pilot=2;%Number of pilot symbols
   

%**********************
% FEC function command
%**********************
strength_len=7;
if strength_len==0
   g=0;
   coding_rate=1;
   filename= 'Uncoded_FSC_10Hz.txt';%File to store the results in
else
   g=g12(strength_len);%generate matrix
   coding_rate=1/2;
   filename= 'Coded_FSC_10Hz.txt';%File to store the results in
   filename1= 'Uncoded_FSC_10Hz.txt';%File to store the results in
end


%*********************************************
% Carriers used for a single wide OFDM channel
%*********************************************
StartCarr=1;
FinCarr=ifftsize;
carriers=[StartCarr:CarrSpacing:FinCarr ];
rep=25;%Number of loop times

SNRMin=0;
SNRMax=30;
SNRInc=5;
perfect=0;%perfect channel option(1: perfect,0:non perfect)
proc_gain=64;

WordSizes=[2];%%BPSK,QPSK,16QAM
NumSizes=length(WordSizes);
Num_Tx=2;%Number of transmiter

headerstr='SNR (samples)';      
BERstr=[];
PhErrstr=[];


%*****************************************************
% Various Delay spread or multipath situation 
% This simulation can simulate time dispersive channel 
% and multipath fading channel
%*****************************************************
Result = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);
Result1 = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);
Result2 = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);
Result3 = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);

%*****************************
% Various modulation schemes.
%*****************************
for l = 1:NumSizes % wordsizeによって変化する
   wordsize = WordSizes(l);
   disp(['wordsize: ' num2str(wordsize)]);
   BERstr = [BERstr ', BER b/Hz: ', int2str(wordsize)];
   PhErrstr = [PhErrstr ', Ph b/Hz: ', int2str(wordsize)];   
   %******************************
   % Repeat various SNR situation
   %******************************
   for k = 1:((SNRMax-SNRMin)/SNRInc)+1
      SNR = (k-1)*SNRInc + SNRMin;
	  disp(['   Eb/N0: ' num2str(SNR)]);
      refData=[]; symwaves=[]; 
      %*******************************
      % Repeat each run several times
      %*******************************
	  for r = 1:rep
        disp(['      Repeat ', num2str(r)]);
        %**********
        %Tx device
        %**********
        Gen_data1=OFDM_data(Num_sym,wordsize,NumCarr);
        Datatx1=transmitter(Gen_data1,Num_sym,wordsize,NumCarr,guardtype,guardtime,...
                            Num_pilot,proc_gain);% [22,64]のデータ（パイロット信号　含）
        
    %**********************************************
	% Make a timewave using IFFT and normalization
	%**********************************************
	%Use example: BaseSignal=(ifft(Datatx'))*sqrt(NumCarr)/sqrt(Num_Tx);%
	BaseSignal1=(ifft(Datatx1'))*sqrt(NumCarr);% [64,22]
    %*******************
	%Add a Guard Period
	%*******************
        if guardtype~=0
           if guardtype == 1%guard period is just a zero transmission
              BaseSignal1=[zeros(guardtime,Numsymb);BaseSignal1];
   		   elseif guardtype == 2
      		  EndSignal = size(BaseSignal1,1);%Find the number of columns in the BaseSignal
      		  BaseSignal1=[BaseSignal1((EndSignal-guardtime+1):EndSignal,:); BaseSignal1]; % [80,22]
           end	
        end
        
        %***************
        %Channel block
        %***************
        %%Rayleigh fading channel
        Multi = zeros(MultiNo,1);
        for w=1:MultiNo
           Multi(w)=10^(0-(((w-1)*ww)/10)); 
        end
        Norm=sum(Multi); NormMulti=Multi./(Norm);%Normalization
        [OutSignal1,Signal_dis1,Fading]=channel_only(BaseSignal1,NormMulti,MultiNo,...
                                        Transrate,Doppler,Delay,SNR,com_delay);
        % channel_only([80,22],[7,1],7,78100*2,10,2,7,0)
        % OutSignal1=[1,80*22]（マルチパス環境では、Multipathfading)
        % Signal_dis1=送信信号[80,22]の標準偏差 ⇒ 1つの値が得られる
        % Fading=[7,22] フェーシング　⇒ 各Pathのシンボルに対してのフェージング
        %%AWGN channel
        TimeSignal1=awgn(OutSignal1,Signal_dis1,MultiNo,Doppler,SNR,Num_sym,...
                         NumCarr,guardtime,Num_pilot,1,perfect,wordsize);% TimeSignal1=[1,22*80]
        clear OutSignal1;%save memory
        clear BaseSignal1;%save memory
        
        %**********
        %Rx device
        %**********
        word=receiver_ams(TimeSignal1,ifftsize,carriers,wordsize,guardtype,...
                          guardtime,Num_sym,Num_pilot,Doppler);%([1,22*80],64,64,1～2,2,16,20,2,10)
                     
        %**********
        %Tx device
        %**********
        if strength_len==0
           Gen_data=OFDM_data(Num_sym,word,NumCarr);
        else
           [Gen_data,Ori_data,int_pattern]=COFDM_data(Num_sym,word,NumCarr,g);
        end
        Datatx=transmitter(Gen_data,Num_sym,word,NumCarr,guardtype,guardtime,...
                           Num_pilot,proc_gain);
                        
    %**********************************************
	% Make a timewave using IFFT and normalization
	%**********************************************
	%Use example: BaseSignal=(ifft(Datatx'))*sqrt(NumCarr)/sqrt(Num_Tx);%
	BaseSignal=(ifft(Datatx'))*sqrt(NumCarr);%        
    %*******************
	%Add a Guard Period
	%*******************
        if guardtype~=0
           if guardtype == 1%guard period is just a zero transmission
              BaseSignal=[zeros(guardtime,Numsymb);BaseSignal];
   		   elseif guardtype == 2
      		  EndSignal = size(BaseSignal,1);%Find the number of columns in the BaseSignal
      		  BaseSignal=[BaseSignal((EndSignal-guardtime+1):EndSignal,:); BaseSignal];
  		   end	
        end
        
        %***************
        %Channel block
        %***************
        [OutSignal,Signal_dis]=close_channel(BaseSignal,NormMulti,MultiNo,Delay,Fading,Doppler);
        %%AWGN channel
        TimeSignal=awgn(OutSignal,Signal_dis,MultiNo,Doppler,SNR,Num_sym,...
                        NumCarr,guardtime,Num_pilot,coding_rate,perfect,word);
        clear OutSignal;%save memory
        clear BaseSignal;%save memory
        
        %**********
        %Rx device
        %**********
        [Datarx,Datarx_hd]=receiver(TimeSignal,ifftsize,carriers,word,guardtype,...
                                    guardtime,Num_sym,Num_pilot,Doppler,proc_gain);

        if strength_len==0
           Summary=error_count(Gen_data,Datarx);
           Result(k,1) = SNR;                  
           %*******************************************
           % Add up the number of errors
           %*******************************************
           Result(k,l+1) = Result(k,l+1)+Summary(1);	
        else
           Summary=error_count(Gen_data,Datarx_hd);
           %[Gen_data(1,1:10,:);Datarx(1,1:10,:)]
           Summary2=error_count2(Ori_data,Datarx,g,int_pattern);
           Result(k,1) = SNR+10*log10(coding_rate);  
           Result1(k,1) = SNR;                  
           Result2(k,1) = SNR;                  
           Result3(k,1) = SNR;                  
           %*******************************************
           % Add up the number of errors
           %*******************************************
           Result(k,l+1) = Result(k,l+1)+Summary(1);	
           Result1(k,l+1) = Result1(k,l+1)+Summary2(1);	
           Result2(k,l+1) = Result2(k,l+1)+word;	
           Result3(k,l+1) = Result3(k,l+1)+Summary2(2);	
        end
     end
   end	
end

%*******************************************************
% Display the result 
%*******************************************************
headerstr = [headerstr, BERstr, PhErrstr];
Result(:,2:NumSizes+1) = Result(:,2:NumSizes+1)/rep;%Average the std dev.
Result1(:,2:NumSizes+1) = Result1(:,2:NumSizes+1)/rep;%Average the std dev.
Result2(:,2:NumSizes+1) = Result2(:,2:NumSizes+1)/rep;%Average the std dev.
Result3(:,2:NumSizes+1) = Result3(:,2:NumSizes+1)/rep;%Average the std dev.
%clear Datatx;clear Datarx;

%*******************************************************
% Display the simulation Time and Processing Speed
%*******************************************************
savefile(filename,Result,headerstr);
if strength_len~=0
   savefile(filename1,Result1,headerstr);
end
disp(['Total Time: ' num2str(toc) 'sec']);
%disp(['Total FLOPS: ' num2str(flops)]);
%disp(['Process Speed : ' num2str(flops/toc) ' flops/sec']);
title('BER ');
plot(SNRMin:SNRInc:SNRMax,(Result2(:,2).*(1-Result3(:,2)))','b<-');%
%semilogy(SNRMin:SNRInc:SNRMax,Result(:,2)','b>-');%,...
%         SNRMin+round(10*log10(coding_rate)):SNRInc:SNRMax...
%               +round(10*log10(coding_rate)),Result1(:,2)','b<-');%
xlabel('Eb/No [dB]');
ylabel('Normalized Throughput [b/Hz]');
%axis([0,30,10^(-5),10^(0)]);
grid;