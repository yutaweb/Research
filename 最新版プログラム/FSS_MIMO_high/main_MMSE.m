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
guardtype=2;% サブキャリアの後方部を前方にコピーする
guardtime=16;% ガードインターバルの長さ
windowtype=2;%　変調方式によって変わる BPSK:1,QPSK:2,etc
CarrSpacing=1;% キャリア間隔 

MultiNo=7;% マルチパスの数
Transrate=20*10^6;
Delay=2;%Number of delay sample
ww=1;%dB degradation ⇒ dB劣化

NumCarr=64;
Num_sym=20;%Number of data symbols    
Num_pilot=2;%Number of pilot symbols
   
%**********************
% FEC function command ⇒ 誤り訂正部
%**********************
strength_len=0;% 拘束長
if strength_len==0
   g=0;
   coding_rate=1;% 符号化率：有用な符号の割合（冗長成分を除く）
   filename= 'result_uncoded_10Hz.txt';%File to store the results in
else
   g=g12(strength_len);% generate matrix
   coding_rate=1/2;% R ⇒ 入力1に対して、出力が2つ得られる。
   filename= 'result_coded_10Hz.txt';%File to store the results in
   filename1= 'result_uncoded_10Hz.txt';%File to store the results in
   filename2= 'result_throughput_10Hz.txt';%File to store the results in
end

%*********************************************
% Carriers used for a single wide OFDM channel
%*********************************************
StartCarr=1;
FinCarr=ifftsize;% FinCarr=64
carriers=[StartCarr:CarrSpacing:FinCarr];% 0～64まで「1」間隔でcarriersにデータが入る [1,64]
rep=1000;%Number of loop times ⇒　任意の回数ループする（シミュレーション回数）

SNRMin=0;
SNRMax=30;
SNRInc=5;
perfect=0;%perfect channel option(1: perfect,0:non perfect)
proc_gain=64;

WordSizes=[2];%%BPSK：1,QPSK：2,16QAM：4
NumSizes=length(WordSizes);% Numsizes=1 ⇒　行列の大きさ
Num_Tx=2;%Number of transmiter ⇒ 送信アンテナの数 （2by2 MIMO）
Num_Rx=2;%Number of receiver ⇒ 受信アンテナの数 （2by2 MIMO）

headerstr='SNR (samples)';      
BERstr=[];
PhErrstr=[];

%*****************************************************
% Various Delay spread or multipath situation 
% This simulation can simulate time dispersive channel 
% and multipath fading channel
%*****************************************************
% ゼロを割り当ててメモリを確保している ⇒ データを記憶する
Result = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)
Result1 = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)
Result2 = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)

%*****************************
% Various modulation schemes.
%*****************************
for l = 1:NumSizes % 1
   wordsize = WordSizes(l);% wordsize = 2
   disp(['wordsize: ' num2str(wordsize)]);% wordsize : 2が表示される
   BERstr = [BERstr ', BER b/Hz: ', int2str(wordsize)];% int2strは整数を文字に変換する
   PhErrstr = [PhErrstr ', Ph b/Hz: ', int2str(wordsize)];   
   %******************************
   % Repeat various SNR situation
   %******************************
   for k = 1:((SNRMax-SNRMin)/SNRInc)+1 % 1～7 
      SNR = (k-1)*SNRInc + SNRMin;% 0～30 ※SNRInc=3
	  disp(['   Eb/N0: ' num2str(SNR)]);% 0～30までを表示する
      refData=[]; symwaves=[]; 
      %*******************************
      % Repeat each run several times
      %*******************************
	  for r = 1:rep % 1～100 ⇒ 各々のEb/Noに対して任意の回数シミュレーションする
        disp(['      Repeat ', num2str(r)]);
        %*********************
        % Tx device ⇒ 送信機
        %*********************
        if strength_len==0
           Gen_data=MIMO_data(Num_sym,wordsize,NumCarr,Num_Tx);% MIMO_data(20,2,64,2)⇒ 最終的に[20,64,2,2]の0or1が返ってくる
        else
           % Original Version
           %[Gen_data,Ori_data,int_pattern]=CMIMO_data(Num_sym,wordsize,NumCarr,g,Num_Tx);% convolutional code(R=1/2,K=7)
           % Gen_data=[20,64,2,2]←0or1,Ori_data=[1,20*64*2-6]←0or1,int_pattern=[1,20*64*2*2]←index
           
           % Fixed Version
           [Gen_data,Ori_data,int_pattern]=CMIMO_data(Num_sym,wordsize,NumCarr,g,Num_Tx);% convolutional code(R=1/2,K=7)
           % Gen_data=[20,64,2,2]←0or1,Ori_data=[1,20*64-6,2]←0or1,int_pattern=[20*64*2,2]←index
        end
        Datatx=transmitter_non_FSS(Gen_data,Num_sym,wordsize,NumCarr,guardtype,guardtime,...
                           Num_pilot,proc_gain,Num_Tx); % transmitter([20,64,2,2],20,2,64,2,16,2,64,2) ⇒ [22,64,2]のデータ（パイロット信号 含）
       
        % [22,64,2]
        % 拡散成分を3次元方向に足し合わせているので、いつもと異なった求め方になっていることに注意する。
        % **********************************************
		%  Make a timewave using IFFT and normalization
		% **********************************************
		% Use example: BaseSignal=(ifft(Datatx'))*sqrt(NumCarr)/sqrt(Num_Tx);%
       
        BaseSignal=zeros(NumCarr,Num_sym+Num_pilot,Num_Tx); % [64,22,2]
        for i=1:Num_Tx % 1～2
            BaseSignal(:,:,i)=(ifft(Datatx(:,:,i)'))*sqrt(NumCarr);% ⇒ [64,22,i]（Datatxで正規化の為に、sqrt(64)で割ったので、ifft部では、sqrt(64)を掛けている）
        end
        % 複素共役転置を取った後にifftをしている。
        % ifft,fftはifft(x,[],1)で列方向のifft(fft)
        %           ifft(x,[],2)で行方向のifft(fft)
        % ※何も指定しない場合は、ifft(x,[],1)と同じ扱いになる。
        % [64,22,2]
        %*******************
		% Add a Guard Period
		%*******************
        if guardtype~=0% 等しくない場合は1を返し、等しい場合は0を返す。
            if guardtype == 1%guard period is just a zero transmission
              BaseSignal=[zeros(guardtime,Num_sym+Num_pilot,Num_Tx);BaseSignal];
            elseif guardtype == 2% ⇒ 64サブキャリアの後方成分を前方に持ってくる
      		  EndSignal = size(BaseSignal,1);% Find the number of columns in the BaseSignal ⇒ 64
      		  BaseSignal=[BaseSignal((EndSignal-guardtime+1):EndSignal,:,:); BaseSignal];% [80,22,2] ⇒ ガードインターバルをつけた
            end	
        end        
        %*****************
        %  Channel block
        %*****************
        %%  Rayleigh fading channel⇒フラットフェージング（時間には依存しない&周波数に対して一様）
        Multi = zeros(MultiNo,1);% [7,1]の0
        for w=1:MultiNo% 1～7
           Multi(w)=10^(0-(((w-1)*ww)/10));% 10^0～10^(-6/10)⇒dBを直した
        end
        Norm=sum(Multi); NormMulti=Multi./(Norm);% Normalization⇒正規化している（Multiの合計値が1になる）
        OutSignal1=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx);
        OutSignal=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx);
        Signal_dis1=zeros(Num_Rx);
        Signal_dis=zeros(Num_Rx);
        Fading1=zeros(MultiNo,Num_sym+Num_pilot,Num_Rx);
        Fading=zeros(MultiNo,Num_sym+Num_pilot,Num_Rx);
        for i=1:Num_Rx % 1～2
            for j=1:Num_Tx % 1～2 Signal_disは後で改善する ⇒ ただデータを保存してるだけなので、Fadingに直接の影響はない。
                [OutSignal1(:,:,j),Signal_dis1(j), Fading1(:,:,j)]=channel_only(BaseSignal(:,:,j),NormMulti,MultiNo,...
                                            Transrate,Doppler,Delay,SNR,com_delay);% channel_only([80,22,2],[7,1],7,1000,4,2,7,0)
                OutSignal(:,:,i)=OutSignal(:,:,i)+OutSignal1(:,:,j);
                Signal_dis(i)=Signal_dis(i)+Signal_dis1(j);
                Fading(:,:,i)=Fading(:,:,i)+Fading1(:,:,j);
            end
        end
        
        % 最初の160成分は一致していてよい。⇒IFFTを掛けているから。
        % OutSignal=[1,80*22,2]（マルチパス環境では、Multipathfading)
        % OutSignal[1,80*22,1]⇒h1*x1+h2'*x2
        % OutSignal[1,80*22,2]⇒h1'*x1+h2*x2
        % ※hはチャネル応答を表す
        % Signal_dis=送信信号[80,22,2]の標準偏差 ⇒ 2つの値が得られる　※若干怪しい..
        % Fading=[7,22,2] フェーシング　⇒ 各Pathのシンボルに対してのフェージング
        %%  AWGN channel
        TimeSignal=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx); % [1,22*80,2]
        noise=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx); % [1,22*80,2]
        for receive=1:Num_Rx
             [TimeSignal(:,:,receive),noise]=awgn(OutSignal(:,:,receive),Signal_dis(receive),MultiNo,Doppler,SNR,Num_sym,...
                           NumCarr,guardtime,Num_pilot,coding_rate,perfect,wordsize);% [[1,22*80,2],1,7,4,7,20,64,16,2,cording_rate,0,1]
        end
        % TimeSignal[1,80*22,1]⇒h1*x1+h2'*x2+noise1
        % TimeSignal[1,80*22,2]⇒h1'*x1+h2*x2+noise2
        % ⇒ 白色雑音を乗せる。[1,22*80,2]
        % 白色雑音：自然界で発生するランダムな雑音
        clear OutSignal;% save memory
        clear BaseSignal;% save memory
  
        %**************
        % Rx device
        %**************
        [Datarx,Datarx_hd]=receiver_non_FSS_MMSE(TimeSignal,ifftsize,carriers,wordsize,guardtype,...
                                    guardtime,Num_sym,Num_pilot,Doppler,proc_gain,Num_Tx,Num_Rx,noise);% [[1,22*80,2],64,[1,64],1,2,4,20,2,4,64,2,2]
        
        % Datarx_hd=復調したデータ [20,64,2,2]
        % Datarx=実数部の生データを符号反転したもの [20,64,2,2]
        
        if strength_len==0
           Summary=error_count(Gen_data,Datarx_hd,Num_Tx);
           Result(k,1) = SNR;                  
           %*******************************************
           % Add up the number of errors
           %*******************************************
           Result(k,l+1) = Result(k,l+1)+Summary(1);	
        else
           Summary=error_count(Gen_data,Datarx_hd,Num_Tx); % Gen_data=[20,64,2,2],Datarx_hd=[20,64,2,2],Num_Tx=2
           Summary2=error_count2(Ori_data,Datarx,g,int_pattern,Num_Tx); % Ori_data=[1,20*64*2-6],Datarx=[20,64,2,2],g=[2,7],int_pattern=[1,20*64*2*2]←index
           Result(k,1) = SNR+10*log10(coding_rate);  % cording_rate=1/2
           Result1(k,1) = SNR;                  
           %*******************************************
           % Add up the number of errors
           %*******************************************
           Result(k,l+1) = Result(k,l+1)+Summary(1);	
           Result1(k,l+1) = Result1(k,l+1)+Summary2(1);	
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
clear Datatx;clear Datarx;

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
semilogy(SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),Result(:,2)','r>-',...
         SNRMin:SNRInc:SNRMax,Result1(:,2)','b<-');
legend('FSS-MIMO(MMSE_MLD)','FSS-MIMO(FEC)');
xlabel('Eb/No [dB]');
ylabel('BER');
axis([0,30,10^(-5),10^(0)]);
grid;
hold on;