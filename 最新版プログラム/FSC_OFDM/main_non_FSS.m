%******************************************************************
% 10/9～ MIMOへの拡張作業をする（遅くとも今週中(土日中)には着手）
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
% flops(0); ⇒　処理の性能を評価するためのものであるので、特に意味はない（現在はtic等で評価）
tic;%Measure the time it takes to run the simulation
rand('state',sum(100*clock));%random probability

%***************************
% Definition for simulation
%***************************
Doppler=10; % 10Hz
com_delay=0;% ??

ifftsize=64;
guardtype=2;% サブキャリアの後方部を前方にコピーする
guardtime=16;% ガードインターバルの長さ
windowtype=2;%　変調方式によって変わる BPSK:1,QPSK:2,etc
CarrSpacing=1;% キャリア間隔 

MultiNo=7;% マルチパスの数
%Transrate=78100*2; % Multipath fadingで用いる
Transrate=20*10^6;
%Translate=1000; 
Delay=2;%Number of delay sample
ww=1;%dB degradation ⇒ dB劣化　(1=>-1dB,2=>-2dB)

NumCarr=64;
Num_sym=20;%Number of data symbols    
Num_pilot=2;%Number of pilot symbols
   
%**********************
% FEC function command ⇒ 誤り訂正部
%**********************
strength_len=7;% K 
if strength_len==0
   g=0;
   coding_rate=1;% 符号化率：有用な符号の割合（冗長成分を除く）
   filename= 'result_uncoded_10Hz.txt';%File to store the results in
else
   g=g12(strength_len);%generate matrix [2,7]
   % g = [1 1 1 1 0 0 1;
   %      1 0 1 1 0 1 1];
   coding_rate=1/2;% R
   filename= 'result_coded_10Hz.txt';%File to store the results in
   filename1= 'result_uncoded_10Hz.txt';%File to store the results in
   filename2= 'result_throughput_10Hz.txt';%File to store the results in
end

%*********************************************
% Carriers used for a single wide OFDM channel
%*********************************************
StartCarr=1;
FinCarr=ifftsize;% FinCarr=64
carriers=[StartCarr:CarrSpacing:FinCarr ];% 0～64まで「1」間隔でcarriersにデータが入る
rep=500;%Number of loop times ⇒　100回ループする（シミュレーション回数）

SNRMin=0;
SNRMax=30;
SNRInc=5;
perfect=0;%perfect channel option(1: perfect,0:non perfect)
%proc_gain=64;% ??
proc_gain=64;

WordSizes=[2];%%BPSK：1,QPSK：2,16QAM：4
NumSizes=length(WordSizes);% Numsizes=1 ⇒　行列の大きさ
Num_Tx=2;%Number of transmiter ⇒ 送信アンテナの数（OFDMは1）ユーザーが増えれば、増やす

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
   wordsize = WordSizes(l);% wordsize = 2 (QPSK)
   disp(wordsize);
   disp(['wordsize: ' num2str(wordsize)]);% wordsize : 2 が表示される
   BERstr = [BERstr ', BER b/Hz: ', int2str(wordsize)];% int2strは整数を文字に変換する
   PhErrstr = [PhErrstr ', Ph b/Hz: ', int2str(wordsize)];   
   %******************************
   % Repeat various SNR situation
   %******************************
   for k = 1:((SNRMax-SNRMin)/SNRInc)+1 % 1～7 
      SNR = (k-1)*SNRInc + SNRMin;% 0～30 ※SNRInc=3
	  disp(['   Eb/N0: ' num2str(SNR)]);% 0～30までを表示する
      refData=[]; symwaves=[]; % ??
      %*******************************
      % Repeat each run several times
      %*******************************
	  for r = 1:rep % 1～100 ⇒ 各々のEb/Noに対して100回シミュレーションする
        disp(['      Repeat ', num2str(r)]);
        %**********
        % Tx device ⇒ 送信機
        %**********
        if strength_len==0 % 拘束長：誤り訂正部
           Gen_data=OFDM_data(Num_sym,wordsize,NumCarr); % OFDM_data(20,1,64)⇒[20,64,1]の0or1が返ってくる
        else
           [Gen_data,Ori_data,int_pattern]=COFDM_data(Num_sym,wordsize,NumCarr,g);% convolutional code(R=1/2,K=7)
           % FECを用いて符号化されたOFDM_dataが入っている。FEC（前方誤り訂正）なので、あらかじめ意図的に符号化したデータを用いている。
           % Gen_data=[20,64,2],int_pattern=[1,20*64*2]⇒データをランダムに並べ替えたときのindex,Ori_data=[1,20*64-6]の0or1
           % Ori_dataをencoderに通してFEC処理を行うとデータが得られる。
           % 上記処理後、データ[20,64,2]の0pr1データ分だけランダムな実数を発生させ、これをsortする。⇒int_pattern,temp_1（この値は使用しない）が求まる。
           % [20,64,2]の0or1データをソートのindexを用いてランダムに並べ替える。
           % 最終的にGen_dataが得られる。
        end
        
        Datatx=transmitter_non_FSS(Gen_data,Num_sym,wordsize,NumCarr,guardtype,guardtime,...
                           Num_pilot,proc_gain); % transmitter([20,64,2],20,2,64,2,16,2,64 ) ⇒ [22,64,2]のデータ（パイロット信号 含）
        % 拡散成分を3次元方向に足し合わせているので、いつもと異なった求め方になっていることに注意する。
        
        %**********************************************
		% Make a timewave using IFFT and normalization
		%**********************************************
		%Use example: BaseSignal=(ifft(Datatx'))*sqrt(NumCarr)/sqrt(Num_Tx);%
		BaseSignal=(ifft(Datatx'))*sqrt(NumCarr);% ⇒ [64,22,2]（Datatxで正規化の為に、sqrt(64)で割ったので、ifft部では、sqrt(64)を掛けている）
        % 複素共役転置を取った後にifftをしている。⇒ 雑音が乗る前なので、BPSKでは複素共役の意味はない
        % ifft,fftはifft(x,[],1)で列方向のifft(fft)
        %           ifft(x,[],2)で行方向のifft(fft)
        % ※何も指定しない場合は、ifft(x,[],1)と同じ扱いになる。
        
        %*******************
		% Add a Guard Period
		%*******************
        if guardtype~=0% 等しくない場合は1を返し、等しい場合は0を返す。
            if guardtype == 1%guard period is just a zero transmission
              BaseSignal=[zeros(guardtime,Numsymb);BaseSignal];
            elseif guardtype == 2% ⇒ 64サブキャリアの後方成分を前方に持ってくる
      		  EndSignal = size(BaseSignal,1);% Find the number of columns in the BaseSignal ⇒ 64
      		  BaseSignal=[BaseSignal((EndSignal-guardtime+1):EndSignal,:); BaseSignal];% [80,22] ⇒ ガードインターバルをつけた
            end	
        end
        
        %***************
        % Channel block
        %***************
        %%  Rayleigh fading channel⇒フラットフェージング（時間には依存しない&周波数に対して一様）
        Multi = zeros(MultiNo,1);% [7,1]の0
        for w=1:MultiNo% 1～7
           Multi(w)=10^(0-(((w-1)*ww)/10));% 10^0～10^(-6/10)⇒dBを直した
        end
        Norm=sum(Multi); NormMulti=Multi./(Norm);%Normalization⇒正規化している（Multiの合計値が1になる）
        [OutSignal,Signal_dis, Fading]=channel_only(BaseSignal,NormMulti,MultiNo,...
                                       Transrate,Doppler,Delay,SNR,com_delay);% channel_only([80,22],[7,1],7,任意の値,4,2,7,0)
      
        % OutSignal=[1,80*22]（マルチパス環境では、Multipathfading)
        % Signal_dis=送信信号[80,22]の標準偏差 ⇒ 1つの値が得られる
        % Fading=[7,22] フェーシング　⇒ 各Pathのシンボルに対してのフェージング
        %%  AWGN channel
        TimeSignal=awgn(OutSignal,Signal_dis,MultiNo,Doppler,SNR,Num_sym,...
                        NumCarr,guardtime,Num_pilot,coding_rate,perfect,wordsize);% [[1,22*80],1,7,4,7,20,64,16,2,1/2,0,1]
        % ⇒ 白色雑音を乗せる。[1,22*80]
        % 白色雑音：自然界で発生するランダムな雑音
        clear OutSignal;%save memory
        clear BaseSignal;%save memory
        
        %**************
        % Rx device
        %**************
        [Datarx,Datarx_hd]=receiver_non_FSS(TimeSignal,ifftsize,carriers,wordsize,guardtype,...
                                    guardtime,Num_sym,Num_pilot,Doppler,proc_gain);% [[1,22*80],64,64,1,2,16,20,2,4,64]
        % Datarx_hd=復調したデータ [20,64,2]
        % Datarx=実数部の生データを符号反転したもの [20,64,2]
                     
        if strength_len==0
           Summary=error_count(Gen_data,Datarx_hd);
           Result(k,1) = SNR;                  
           %*******************************************
           % Add up the number of errors
           %*******************************************
           Result(k,l+1) = Result(k,l+1)+Summary(1);	
        else
           Summary=error_count(Gen_data,Datarx_hd); % 元の生成データと比較して誤っているデータ数を検出する。⇒ BERが返却される。
           Summary2=error_count2(Ori_data,Datarx,g,int_pattern); 
           % Ori_data=[1,20*64-6]の0or1,Datarx=[20,64,2] ⇒ 0or1に復調していない実数成分（符号反転している）,g=[2,7],int_pattern=[1,20*64*2]のindex
           Result(k,1) = SNR+10*log10(coding_rate); % [20,64,2]のデータを評価している（いつものOFDM） ⇒ Result1()と性能を比較するために、調整している。
           Result1(k,1) = SNR;% 20*64-6のデータを評価している。
           % 軟判定ビタビアルゴリズム、FEC、インターリーブを利用している。
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
legend('FSS-OFDM(2)','FSS-OFDM(FEC)');
xlabel('Eb/No [dB]');
ylabel('BER');
axis([0,30,10^(-5),10^(0)]);
grid;
hold on;
