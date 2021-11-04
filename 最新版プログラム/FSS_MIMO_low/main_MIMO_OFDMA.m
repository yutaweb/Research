% 全体的な課題：MIMO-OFDMAやそもそものOFDMAに関する知識をもっとつける
% チャネルのフィードバックや移動通信環境下での評価をもっと行う
% FSSで受信電力が等しくなっているかを検証する。
% スループットの評価をする
% 高次元変調に対応させる
% SNR(高)：性能が良くなる
% SNR(低)：性能が悪くなる
% ※とりあえずは、16ビットで行う(NG)→なぜ？
% ※4G：64ビット（64QAM）
% ※5G：256ビット（256QAM）
% 最適化の部分の見直し(OK)
% 直交性補償の見直し
% 誤り訂正符号を導入した際のシステムの評価(OK)
% アンテナの本数を変えて評価する
% 電力の見直し(OK)
% 関数名の見直し(OK)
% OFDMAでは、完全に分離できている前提ではなく、SICなどで分離してユーザー間干渉を考慮した上で、評価する必要がある。
% もしくは、ユーザ毎に異なる拡散符号を割り当てる手法を考える。

%% simuration start
clear all;
tic;
rand('state',sum(100*clock));

%% definition of parametor
Doppler=10;
com_delay=0;
ifftsize=64;
guardtype=2;
guardtime=16;
windowtype=2;% BPSK:1,QPSK:2,etc
CarrSpacing=1;
MultiNo=7;% multipath-fading
Transrate=20*10^6;
Delay=2;%Number of delay sample
ww=1;%dB degradation

NumCarr=64;
Num_sym=20;%Number of data symbols    
Num_pilot=2;%Number of pilot symbols
   
%% FEC Function Command

strength_len=0;
if strength_len==0
   g=0;
   coding_rate=1;
   filename= 'result_uncoded_10Hz.txt';%File to store the results in
else
   g=g12(strength_len);% generate matrix
   coding_rate=1/2;
   filename= 'result_coded_10Hz.txt';%File to store the results in
   filename1= 'result_uncoded_10Hz.txt';%File to store the results in
   filename2= 'result_throughput_10Hz.txt';%File to store the results in
end

%% Carriers used for a single wide OFDM channel

StartCarr=1;
FinCarr=ifftsize;% FinCarr =64
carriers=[StartCarr:CarrSpacing:FinCarr];% [1,64]
rep=100;%Number of loop times 

SNRMin=0;
SNRMax=30;
SNRInc=5;
perfect=0;%perfect channel option(1: perfect,0:non perfect)
proc_gain=16; %OFDMA(subcarrier num per 1 block)
WordSizes=2;%%BPSK:1,QPSK:2,16QAM:4
NumSizes=length(WordSizes);% Numsizes=1
Num_Tx=2;%Number of transmiter 
Num_Rx=2;%Number of receiver
Num_User=4;%Number of user

rng('default');
rng(10);
channel_rand=rand([Num_Tx,Num_Rx,Num_User,2,MultiNo,(SNRMax-SNRMin)/SNRInc+1,rep]); 
% [2,2,4,7,7,10]

headerstr='SNR (samples)';      
BERstr=[];
PhErrstr=[];

Result = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)
Result1 = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)
Result2 = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)

%% Various modulation schemes.

for l = 1:NumSizes % 1
   wordsize = WordSizes(l);% wordsize = 2
   disp(['wordsize: ' num2str(wordsize)]);% wordsize : 2
   BERstr = [BERstr ', BER b/Hz: ', int2str(wordsize)];% int2str
   PhErrstr = [PhErrstr ', Ph b/Hz: ', int2str(wordsize)];   
  
   %% Repeat various SNR situation
  
   for k = 1:((SNRMax-SNRMin)/SNRInc)+1 % 1 - 7 
      SNR = (k-1)*SNRInc + SNRMin;% 0 - 30 SNRInc=3
	  disp(['   Eb/N0: ' num2str(SNR)]);% 0 - 30
      refData=[]; symwaves=[]; 
      
      %% Repeat each run several times
      for r = 1:rep % 1 - 100
        disp(['      Repeat ', num2str(r)]);
        %% Tx device
        if strength_len==0
           Gen_data=MIMO_data(Num_sym,wordsize,NumCarr,Num_Tx);% MIMO_data(20,2,64,2)
        else
           [Gen_data,Ori_data,int_pattern]=CMIMO_data(Num_sym,wordsize,NumCarr,g,Num_Tx);
           % convolutional code(R=1/2,K=7)
           % Gen_data=[20,64,2,2] 0or1,Ori_data=[1,20*64*2-6] 0or1,int_pattern=[1,20*64*2*2] index
        end
        
        %% estimate channel
        [channel_ranking,user_index]=main_virtual(SNR,Doppler,com_delay,ifftsize,guardtype,guardtime,...
            MultiNo,Transrate,Delay,ww,Num_sym,NumCarr,Num_pilot,strength_len,carriers,proc_gain,...
            wordsize,NumSizes,k,Num_Tx,Num_Rx,Num_User,channel_rand,r);
       
        %% transmitter
         Datatx=transmitter_non_FSS_MIMO_OFDMA(Gen_data,Num_sym,wordsize,NumCarr,guardtype,guardtime,...
                          Num_pilot,proc_gain,Num_Tx); % transmitter([20,64,2,2],20,2,64,2,16,2,64,4) [22,64,2]
        % [22,64,2]
        % [22,1:16,2]
        % [22,17:32,2]
        % [22,33:48,2]
        % [22,49:64,2]
       %%  Make a timewave using IFFT and normalization
        BaseSignal=zeros(NumCarr,Num_sym+Num_pilot,Num_Tx); % [64,22,2]
        for i=1:Num_Tx % 1 - 2
            BaseSignal(:,:,i)=(ifft(Datatx(:,:,i)'))*sqrt(NumCarr);% [64,22,i]
        end
        
        %% Add a Guard Period
		if guardtype~=0
            if guardtype == 1%guard period is just a zero transmission
              BaseSignal=[zeros(guardtime,Num_sym+Num_pilot,Num_Tx);BaseSignal];
            elseif guardtype == 2
      		  EndSignal = size(BaseSignal,1);% Find the number of columns in the BaseSignal 64
      		  BaseSignal=[BaseSignal((EndSignal-guardtime+1):EndSignal,:,:); BaseSignal];% [80,22,2]
            end	
        end        
        %%  Rayleigh fading channel
        Multi = zeros(MultiNo,1);% [7,1]
        for w=1:MultiNo% 1 - 7
           Multi(w)=10^(0-(((w-1)*ww)/10));% 10^0 - 10^(-6/10)
        end
        Norm=sum(Multi); NormMulti=Multi./(Norm);% Normalization
        OutSignal1=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        OutSignal=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        Signal_dis1=zeros(Num_Rx,Num_User); % [2,4]
        Signal_dis=zeros(Num_Rx,Num_User); % [2,4]
        Fading1=zeros(MultiNo,Num_sym+Num_pilot,Num_Rx,Num_User); % [7,22,2,4]
        Fading=zeros(MultiNo,Num_sym+Num_pilot,Num_Rx,Num_User); % [7,22,2,4]
        for u=1:Num_User % 1 - 4
            for i=1:Num_Rx % 1 - 2
                for j=1:Num_Tx % 1 - 2
                    [OutSignal1(:,:,j,u),Signal_dis1(j,u), Fading1(:,:,j,u)]=channel_only_MIMO_OFDMA(BaseSignal(:,:,j),NormMulti,MultiNo,...
                    Transrate,Doppler,Delay,SNR,com_delay,channel_rand(j,i,u,:,:,k,r));% channel_only([80,22,2],[7,1],7,1000,4,2,7,0)
                    OutSignal(:,:,i,u)=OutSignal(:,:,i,u)+OutSignal1(:,:,j,u); % [1,1760,2,4]
                    Signal_dis(i,u)=Signal_dis(i,u)+Signal_dis1(j,u); % [2,4]
                    Fading(:,:,i,u)=Fading(:,:,i,u)+Fading1(:,:,j,u); % [7,22,2,4]
                end
            end    
        end
        % OutSignal=[1,80*22,2]
        % OutSignal[1,80*22,1] h1*x1+h2'*x2
        % OutSignal[1,80*22,2] h1'*x1+h2*x2
        % Signal_dis=[80,22,2]
        % Fading=[7,22,2]
        %%  AWGN channel
        TimeSignal=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        noise=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        for u=1:Num_User
            for receive=1:Num_Rx
                [TimeSignal(:,:,receive,u),noise(1,:,receive,u)]=awgn(OutSignal(:,:,receive,u),Signal_dis(receive,u),MultiNo,Doppler,SNR,Num_sym,...
                    NumCarr,guardtime,Num_pilot,coding_rate,perfect,wordsize);% [[1,22*80,2,4],1,7,4,7,20,64,16,2,cording_rate,0,1]
            end
        end

        % TimeSignal[1,80*22,1]h1*x1+h2'*x2+noise1
        % TimeSignal[1,80*22,2]h1'*x1+h2*x2+noise2
        clear OutSignal;% save memory
        clear BaseSignal;% save memory
  
         %% Rx device 
        [Datarx,Datarx_hd]=receiver_non_FSS_MIMO_OFDMA(TimeSignal,ifftsize,carriers,wordsize,guardtype,...
                                    guardtime,Num_sym,Num_pilot,Doppler,proc_gain,Num_Tx,Num_Rx,Num_User,noise,channel_ranking,user_index);% [[1,22*80,2,4],64,[1,64],1,2,4,20,2,4,64,2,2,[1,22*80,2,4]]
        
        if strength_len==0
           Summary=error_count(Gen_data,Datarx_hd,Num_Tx);
           Result(k,1) = SNR;
            %% Add up the number of errors
           Result(k,l+1) = Result(k,l+1)+Summary(1);
        else
           Summary=error_count(Gen_data,Datarx_hd,Num_Tx); % Gen_data=[20,64,2,2],Datarx_hd=[20,64,2,2],Num_Tx=2
           Summary2=error_count2(Ori_data,Datarx,g,int_pattern,Num_Tx); % Ori_data=[1,20*64*2-6],Datarx=[20,64,2,2],g=[2,7],int_pattern=[1,20*64*2*2] index
           Result(k,1) = SNR+10*log10(coding_rate);  % cording_rate=1/2
           Result1(k,1) = SNR;                  
            %% Add up the number of errors
           Result(k,l+1) = Result(k,l+1)+Summary(1);	
           Result1(k,l+1) = Result1(k,l+1)+Summary2(1);	
        end
     end
   end	
end

%% Display the result 
headerstr = [headerstr, BERstr, PhErrstr];
Result(:,2:NumSizes+1) = Result(:,2:NumSizes+1)/rep;%Average the std dev.
Result1(:,2:NumSizes+1) = Result1(:,2:NumSizes+1)/rep;%Average the std dev.
clear Datatx;clear Datarx;

%% Display the simulation Time and Processing Speed
savefile(filename,Result,headerstr);
if strength_len~=0
   savefile(filename1,Result1,headerstr);
end
disp(['Total Time: ' num2str(toc) 'sec']);
%disp(['Total FLOPS: ' num2str(flops)]);
%disp(['Process Speed : ' num2str(flops/toc) ' flops/sec']);
semilogy(SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),Result(:,2)','r>-',...
         SNRMin:SNRInc:SNRMax,Result1(:,2)','b<-');
legend('FSS-MIMO-OFDMA','FSS-MIMO-OFDMA(FEC)');
xlabel('Eb/No [dB]');
ylabel('BER');
axis([0,30,10^(-5),10^(0)]);
grid;
hold on;