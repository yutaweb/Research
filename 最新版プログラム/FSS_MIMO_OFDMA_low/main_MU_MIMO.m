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

%*****************
% 課�?
% �?2�?2MIMOを用�?る�?�合�?�ユーザ数Num_Userとproc_gainを変更する�?けで、良�?プログラ�?にする
% ②そ�?�他機�?�をまとめる
% ③SISOと4�?4MIMOにつ�?ては時間があれ�?�、検討す�?
%*****************
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
guardtype=2;% サブキャリアの後方部を前方にコピ�?�する
guardtime=16;% ガードインターバルの長�?
windowtype=2;%�?変調方式によって変わ�? BPSK:1,QPSK:2,etc
CarrSpacing=1;% キャリア間隔 

MultiNo=7;% マルチパスの数
%Transrate=1000;% Multipath fadingで用�?�?
%Transrate=78100*2;
Transrate=20*10^6;
Delay=2;%Number of delay sample
ww=1;%dB degradation �? dB劣�?

NumCarr=64;
Num_sym=20;%Number of data symbols    
Num_pilot=2;%Number of pilot symbols
   
%**********************
% FEC function command �? 誤り訂正部
%**********************
strength_len=0;% 拘束長
if strength_len==0
   g=0;
   coding_rate=1;% 符号化率?��有用な符号の割合（�?�長成�?を除く�?
   filename= 'result_uncoded_10Hz.txt';%File to store the results in
else
   g=g12(strength_len);% generate matrix
   coding_rate=1/2;% R �? 入�?1に対して、�?�力が2つ得られる�?
   filename= 'result_coded_10Hz.txt';%File to store the results in
   filename1= 'result_uncoded_10Hz.txt';%File to store the results in
   filename2= 'result_throughput_10Hz.txt';%File to store the results in
end

%*********************************************
% Carriers used for a single wide OFDM channel
%*********************************************
StartCarr=1;
FinCarr=ifftsize;% FinCarr=64
carriers=[StartCarr:CarrSpacing:FinCarr];% 0?�?64まで�?1」間隔でcarriersに�?ータが�?��? [1,64]
rep=1000;%Number of loop times ⇒�??任意�?�回数ループする（シミュレーション回数?�?

SNRMin=0;
SNRMax=30;
SNRInc=5;
perfect=0;%perfect channel option(1: perfect,0:non perfect)
proc_gain=8; %OFDMAによって�?4ユーザに割り当てる�?�で16�?16のマトリ�?クス
WordSizes=[2];%%BPSK?�?1,QPSK?�?2,16QAM?�?4
NumSizes=length(WordSizes);% Numsizes=1 ⇒�??行�?��?�大きさ
Num_Tx=2;%Number of transmiter �? 送信アン�?ナ�?�数 
Num_Rx=2;%Number of receiver �? 受信アン�?ナ�?�数
Num_User=8;%Number of user �? ユーザの数

rng('default');
rng(10); % ラン�?�?の値を変更した�?とき�?�、この数字を変更する
% 条件を統�?させるために先にすべてのラン�?�?値を生成す�?
channel_rand=rand([Num_Tx,Num_Rx,Num_User,2,MultiNo,(SNRMax-SNRMin)/SNRInc+1,rep]); 
% [2,2,4,7,7,rep] �? [送信アン�?ナ数,受信アン�?ナ数,ユーザー数,channel_only�?での使用数,チャネル伝播路数,マルチパス数,SNR反復数,全反復数]

headerstr='SNR (samples)';      
BERstr=[];
PhErrstr=[];

%*****************************************************
% Various Delay spread or multipath situation 
% This simulation can simulate time dispersive channel 
% and multipath fading channel
%*****************************************************
% ゼロを割り当ててメモリを確保して�?�? �? �?ータを記�?�する
Result = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)
Result1 = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)
Result2 = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)

%*****************************
% Various modulation schemes.
%*****************************
for l = 1:NumSizes % 1
   wordsize = WordSizes(l);% wordsize = 2
   disp(['wordsize: ' num2str(wordsize)]);% wordsize : 2が表示され�?
   BERstr = [BERstr ', BER b/Hz: ', int2str(wordsize)];% int2strは整数を文字に変換する
   PhErrstr = [PhErrstr ', Ph b/Hz: ', int2str(wordsize)];   
   %******************************
   % Repeat various SNR situation
   %******************************
   for k = 1:((SNRMax-SNRMin)/SNRInc)+1 % 1?�?7 
      SNR = (k-1)*SNRInc + SNRMin;% 0?�?30 ※SNRInc=3
	  disp(['   Eb/N0: ' num2str(SNR)]);% 0?�?30までを表示する
      refData=[]; symwaves=[]; 
      %*******************************
      % Repeat each run several times
      %*******************************
	  for r = 1:rep % 1?�?100 �? �?�?のEb/Noに対して任意�?�回数シミュレーションする
        disp(['      Repeat ', num2str(r)]);
        %*********************
        % Tx device �? 送信�?
        % OFDMAの場合�?�各ストリー�?からは、同じストリー�?を�?�信し�?�受信側で自�?に割り当てられて�?るデータを取り�?�す�??
        % �?チャネル数�?1として、各送信アン�?ナからストリー�?を�?�信する�?
        % ②チャネル数を増やし�?��?�差�?新アン�?ナからストリー�?を�?�信する�?
        %*********************
        if strength_len==0
           Gen_data=MIMO_data(Num_sym,wordsize,NumCarr,Num_Tx);% MIMO_data(20,2,64,2)�? �?終的に[20,64,2,2]の0or1が返ってくる
        else
           [Gen_data,Ori_data,int_pattern]=CMIMO_data(Num_sym,wordsize,NumCarr,g,Num_Tx);% convolutional code(R=1/2,K=7)
           % Gen_data=[20,64,2,2]�?0or1,Ori_data=[1,20*64*2-6]�?0or1,int_pattern=[1,20*64*2*2]←index
           % FEC及�?�インタリーブを掛ける際は、MU_MIMOの場合�?��?�ストリー�?に同じも�?�を掛けてよいも�?�と思われる�?
        end
        % 先にチャネル推定をして、割り当てるサブキャリアを決める�?
        
        Datatx=transmitter_non_FSS_MU_MIMO(Gen_data,Num_sym,wordsize,NumCarr,guardtype,guardtime,...
                           Num_pilot,proc_gain,Num_Tx); % transmitter([20,64,2,2],20,2,64,2,16,2,64,4) �? [22,64,2]の�?ータ?��パイロ�?ト信号 含?�?
        % �?ストリー�?において64サブキャリアをユーザ数�?�?け�?割して、各ユーザの�?ータに対して拡散符号を割り当て�?
        % [22,64,2]
        % 4ユーザの場�?
        % [22,1 ?�?16,2] �? ユーザ1に割り当て�?
        % [22,17?�?32,2] �? ユーザ2に割り当て�?
        % [22,33?�?48,2] �? ユーザ3に割り当て�?
        % [22,49?�?64,2] �? ユーザ4に割り当て�?
        % 8ユーザの場合�?�サブキャリア�?8�?割�?16ユーザの場合�?�サブキャリア�?16�?割する�?
        % ※ここで、注意しなければならな�?のは、IEEEの規�?�では、サブキャリア間にヌルを設ける�?要があるので、実際に使えるのは64サブキャリアの場合�?��?52�?
        % 拡散成�?�?3次�?方向に足し合わせて�?る�?�で、従来と異なった求め方になって�?ることに注意する�??
        % **********************************************
		%  Make a timewave using IFFT and normalization
		% **********************************************
		% Use example: BaseSignal=(ifft(Datatx'))*sqrt(NumCarr)/sqrt(Num_Tx);%
       
        % 全サブキャリアに�?気にIFFTを掛ける場�? M
        BaseSignal=zeros(NumCarr,Num_sym+Num_pilot,Num_Tx); % [64,22,2]
        for i=1:Num_Tx % 1?�?2
            BaseSignal(:,:,i)=(ifft(Datatx(:,:,i)'))*sqrt(NumCarr);% �? [64,22,i]?�?Datatxで正規化の為に、sqrt(64)で割った�?�で、ifft部では、sqrt(64)を掛けて�?る�?
        end
        
     
        % �?�?共役転置を取った後にifftをして�?る�??
        % ifft,fftはifft(x,[],1)で列方向�?�ifft(fft)
        %           ifft(x,[],2)で行方向�?�ifft(fft)
        % ※何も�?定しな�?場合�?�、ifft(x,[],1)と同じ扱�?になる�??
        % [64,22,2]
        
        %*******************
		% Add a Guard Period
		%*******************
        if guardtype~=0% 等しくな�?場合�?�1を返し、等し�?場合�?�0を返す�?
            if guardtype == 1%guard period is just a zero transmissionリー�?
              BaseSignal=[zeros(guardtime,Num_sym+Num_pilot,Num_Tx);BaseSignal];
            elseif guardtype == 2% �? 64サブキャリアの後方成�?を前方に持ってくる
      		  EndSignal = size(BaseSignal,1);% Find the number of columns in the BaseSignal �? 64
      		  BaseSignal=[BaseSignal((EndSignal-guardtime+1):EndSignal,:,:); BaseSignal];% [80,22,2] �? ガードインターバルをつけた
            end	
        end        
        %*****************
        %  Channel block �? ユーザの数を�??慮する�?要があることに注意す�?
        %*****************
        %%  Rayleigh fading channel⇒フラ�?トフェージング?��時間には依存しな�?&周波数に対して�?様�?
        Multi = zeros(MultiNo,1);% [7,1]の0
        for w=1:MultiNo% 1?�?7
           Multi(w)=10^(0-(((w-1)*ww)/10));% 10^0?�?10^(-6/10)⇒dBを直した
        end
        Norm=sum(Multi); NormMulti=Multi./(Norm);% Normalization⇒正規化して�?る�?Multiの合計�?��?1になる�?
        % 以下�?�ユーザ数�?�?け作る
        OutSignal1=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        OutSignal=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        Signal_dis1=zeros(Num_Rx,Num_User); % [2,4]
        Signal_dis=zeros(Num_Rx,Num_User); % [2,4]
        Fading1=zeros(MultiNo,Num_sym+Num_pilot,Num_Rx,Num_User); % [7,22,2,4]
        Fading=zeros(MultiNo,Num_sym+Num_pilot,Num_Rx,Num_User); % [7,22,2,4]
        for u=1:Num_User % 1?�?4 �? �?ユーザごとに異なる伝搬環�?及�?�受信アン�?ナノイズを持って�?ることに起�?する
            for i=1:Num_Rx % 1?�?2
                for j=1:Num_Tx % 1?�?2 Signal_disは後で改�?する �? ただ�?ータを保存してるだけなので、Fadingに直接の影響はな�?�?
                    [OutSignal1(:,:,j,u),Signal_dis1(j,u), Fading1(:,:,j,u)]=channel_only_AM(BaseSignal(:,:,j),NormMulti,MultiNo,...
                        Transrate,Doppler,Delay,SNR,com_delay,channel_rand(j,i,u,:,:,k,r));% channel_only([80,22,2],[7,1],7,1000,4,2,7,0)
                    OutSignal(:,:,i,u)=OutSignal(:,:,i,u)+OutSignal1(:,:,j,u); % [1,1760,2,4]
                    Signal_dis(i,u)=Signal_dis(i,u)+Signal_dis1(j,u); % [2,4]
                    Fading(:,:,i,u)=Fading(:,:,i,u)+Fading1(:,:,j,u); % [7,22,2,4]
                end
            end    
        end
        % �?初�?�160成�?は�?致して�?てよい。��IFFTを掛けて�?るから�??
        % OutSignal=[1,80*22,2]?���?�ルチパス環�?では、Multipathfading)
        % OutSignal[1,80*22,1]⇒h1*x1+h2'*x2
        % OutSignal[1,80*22,2]⇒h1'*x1+h2*x2
        % ※hはチャネル応答を表�?
        % Signal_dis=送信信号[80,22,2]の標準偏差 �? 2つの値が得られる�?※若干怪しい..
        % Fading=[7,22,2] フェーシング�?�? 各Pathのシンボルに対してのフェージング
        %%  AWGN channel
        TimeSignal=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        noise=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        % 全ユーザの受信アン�?ナ数の合計�?�8本なので、それぞれに異なるノイズを加える�?要がある
        for u=1:Num_User
            for receive=1:Num_Rx
                [TimeSignal(:,:,receive,u),noise(1,:,receive,u)]=awgn(OutSignal(:,:,receive,u),Signal_dis(receive,u),MultiNo,Doppler,SNR,Num_sym,...
                    NumCarr,guardtime,Num_pilot,coding_rate,perfect,wordsize);% [[1,22*80,2,4],1,7,4,7,20,64,16,2,cording_rate,0,1]
            end
        end

        % TimeSignal[1,80*22,1]⇒h1*x1+h2'*x2+noise1
        % TimeSignal[1,80*22,2]⇒h1'*x1+h2*x2+noise2
        % �? 白色雑音を乗せる�??[1,22*80,2]
        % 白色雑音?���?�然界で発生するラン�?�?な雑音
        clear OutSignal;% save memory
        clear BaseSignal;% save memory
  
        %**************
        % Rx device 
        %**************
        [Datarx,Datarx_hd]=receiver_non_FSS_MU_MIMO(TimeSignal,ifftsize,carriers,wordsize,guardtype,...
                                    guardtime,Num_sym,Num_pilot,Doppler,proc_gain,Num_Tx,Num_Rx,Num_User,noise);% [[1,22*80,2,4],64,[1,64],1,2,4,20,2,4,64,2,2,[1,22*80,2,4]]
        
        % Datarx_hd=復調した�?ータ [20,64,2,2]
        % Datarx=実数部の生データを符号反転したも�?� [20,64,2,2]
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