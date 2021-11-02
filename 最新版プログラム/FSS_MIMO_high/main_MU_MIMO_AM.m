%****************************
%　作成者：遠山裕太
%  修正日：2020/12/02
%****************************

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
MultiNo=7;% マルチパス数
Transrate=20*10^6;
Delay=2;%Number of delay sample
ww=1;%dB degradation ⇒ dB劣化

NumCarr=64;
Num_sym=20;%Number of data symbols    
Num_pilot=2;%Number of pilot symbols
   
%% FEC Function Command

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

%% Carriers used for a single wide OFDM channel

StartCarr=1;
FinCarr=ifftsize;% FinCarr=64
carriers=[StartCarr:CarrSpacing:FinCarr];% 0～64まで「1」間隔でcarriersにデータが入る [1,64]
rep=1000;%Number of loop times ⇒　任意の回数ループする（シミュレーション回数）

SNRMin=0;
SNRMax=30;
SNRInc=5;
perfect=0;%perfect channel option(1: perfect,0:non perfect)
proc_gain=16; %OFDMAによって、4ユーザに割り当てるので16×16のマトリックス
WordSizes=[2];%%BPSK：1,QPSK：2,16QAM：4
NumSizes=length(WordSizes);% Numsizes=1 ⇒　行列の大きさ
Num_Tx=2;%Number of transmiter ⇒ 送信アンテナの数 
Num_Rx=2;%Number of receiver ⇒ 受信アンテナの数
Num_User=4;%Number of user ⇒ ユーザの数

% 条件を統一させるために先にすべてのランダム値を生成する
rng('default');
rng(10); % ランダムの値を変更したい時は、この数字を変更する
channel_rand=rand([Num_Tx,Num_Rx,Num_User,2,MultiNo,(SNRMax-SNRMin)/SNRInc+1,rep]); 
% [2,2,4,7,7,10] → [送信アンテナ数,受信アンテナ数,ユーザー数,channel_only内での使用数,チャネル伝播路数,マルチパス数,SNR反復数,全反復数]

headerstr='SNR (samples)';      
BERstr=[];
PhErrstr=[];

Result = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)
Result1 = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)
Result2 = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)

%% Various modulation schemes.

for l = 1:NumSizes % 1
   wordsize = WordSizes(l);% wordsize = 2
   disp(['wordsize: ' num2str(wordsize)]);% wordsize : 2が表示される
   BERstr = [BERstr ', BER b/Hz: ', int2str(wordsize)];% int2strは整数を文字に変換する
   PhErrstr = [PhErrstr ', Ph b/Hz: ', int2str(wordsize)];   
  
   %% Repeat various SNR situation
  
   for k = 1:((SNRMax-SNRMin)/SNRInc)+1 % 1～7 
      SNR = (k-1)*SNRInc + SNRMin;% 0～30 ※SNRInc=3
	  disp(['   Eb/N0: ' num2str(SNR)]);% 0～30までを表示する
      refData=[]; symwaves=[]; 
      
      %% Repeat each run several times
      for r = 1:rep % 1～100 ⇒ 各々のEb/Noに対して任意の回数シミュレーションする
        disp(['      Repeat ', num2str(r)]);
        %% Tx device
        if strength_len==0
           Gen_data=MIMO_data(Num_sym,wordsize,NumCarr,Num_Tx);% MIMO_data(20,2,64,2)⇒ 最終的に[20,64,2,2]の0or1が返ってくる
        else
           [Gen_data,Ori_data,int_pattern]=CMIMO_data(Num_sym,wordsize,NumCarr,g,Num_Tx);% convolutional code(R=1/2,K=7)
           % Gen_data=[20,64,2,2]←0or1,Ori_data=[1,20*64*2-6]←0or1,int_pattern=[1,20*64*2*2]←index
           % FEC及びインタリーブを掛ける際は、MU_MIMOの場合、全ストリームに同じものを掛けてよいものと思われる。
        end
        
        %% estimate channel
        [channel_ranking,user_index]=main_virtual(SNR,Doppler,com_delay,ifftsize,guardtype,guardtime,...
            MultiNo,Transrate,Delay,ww,Num_sym,NumCarr,Num_pilot,strength_len,carriers,proc_gain,...
            WordSizes,NumSizes,k,Num_Tx,Num_Rx,Num_User,channel_rand,r);
       
        %% transmitter
         Datatx=transmitter_non_FSS_MU_MIMO(Gen_data,Num_sym,wordsize,NumCarr,guardtype,guardtime,...
                          Num_pilot,proc_gain,Num_Tx,block_num_per_user); % transmitter([20,64,2,2],20,2,64,2,16,2,64,4) ⇒ [22,64,2]のデータ（パイロット信号 含）
        
        % 各ストリームにおいて64サブキャリアをユーザ数分だけ分割して、各ユーザのデータに対して拡散符号を割り当てる
        % [22,64,2]
        % 4ユーザの場合
        % [22,1 ～16,2] ⇒ ユーザ1に割り当てる
        % [22,17～32,2] ⇒ ユーザ2に割り当てる
        % [22,33～48,2] ⇒ ユーザ3に割り当てる
        % [22,49～64,2] ⇒ ユーザ4に割り当てる
        % 8ユーザの場合はサブキャリアを8分割、16ユーザの場合はサブキャリアを16分割する。
        % ※ここで、注意しなければならないのは、IEEEの規格では、サブキャリア間にヌルを設ける必要があるので、実際に使えるのは64サブキャリアの場合は、52個
        % 拡散成分を3次元方向に足し合わせているので、従来と異なった求め方になっていることに注意する。
        %%  Make a timewave using IFFT and normalization
        % 全サブキャリアに一気にIFFTを掛ける場合 M
        BaseSignal=zeros(NumCarr,Num_sym+Num_pilot,Num_Tx); % [64,22,2]
        for i=1:Num_Tx % 1～2
            BaseSignal(:,:,i)=(ifft(Datatx(:,:,i)'))*sqrt(NumCarr);% ⇒ [64,22,i]（Datatxで正規化の為に、sqrt(64)で割ったので、ifft部では、sqrt(64)を掛けている）
        end
         
        % 各ユーザに割り当てたサブキャリアブロックごとにIFFTを掛ける場合
%         proc_gain_block=NumCarr/proc_gain; % 4
%         for fscb=1:proc_gain_block% 1～4
%             for i=1:Num_Tx % 1～2
%                 BaseSignal((fscb-1)*proc_gain+1:fscb*proc_gain,:,i)=(ifft(Datatx(:,(fscb-1)*proc_gain+1:fscb*proc_gain,i)'))*sqrt(proc_gain);
%             end
%         end
     
        % 複素共役転置を取った後にifftをしている。
        % ifft,fftはifft(x,[],1)で列方向のifft(fft)
        %           ifft(x,[],2)で行方向のifft(fft)
        % ※何も指定しない場合は、ifft(x,[],1)と同じ扱いになる。
        % [64,22,2]
        
        %% Add a Guard Period
		if guardtype~=0% 等しくない場合は1を返し、等しい場合は0を返す。
            if guardtype == 1%guard period is just a zero transmission
              BaseSignal=[zeros(guardtime,Num_sym+Num_pilot,Num_Tx);BaseSignal];
            elseif guardtype == 2% ⇒ 64サブキャリアの後方成分を前方に持ってくる
      		  EndSignal = size(BaseSignal,1);% Find the number of columns in the BaseSignal ⇒ 64
      		  BaseSignal=[BaseSignal((EndSignal-guardtime+1):EndSignal,:,:); BaseSignal];% [80,22,2] ⇒ ガードインターバルをつけた
            end	
        end        
        %%  Rayleigh fading channel⇒フラットフェージング（時間には依存しない&周波数に対して一様）
        Multi = zeros(MultiNo,1);% [7,1]の0
        for w=1:MultiNo% 1～7
           Multi(w)=10^(0-(((w-1)*ww)/10));% 10^0～10^(-6/10)⇒dBを直した
        end
        Norm=sum(Multi); NormMulti=Multi./(Norm);% Normalization⇒正規化している（Multiの合計値が1になる）
        % 以下、ユーザ数分だけ作る
        OutSignal1=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        OutSignal=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        Signal_dis1=zeros(Num_Rx,Num_User); % [2,4]
        Signal_dis=zeros(Num_Rx,Num_User); % [2,4]
        Fading1=zeros(MultiNo,Num_sym+Num_pilot,Num_Rx,Num_User); % [7,22,2,4]
        Fading=zeros(MultiNo,Num_sym+Num_pilot,Num_Rx,Num_User); % [7,22,2,4]
        for u=1:Num_User % 1～4 ⇒ 各ユーザごとに異なる伝搬環境及び受信アンテナノイズを持っていることに起因する
            for i=1:Num_Rx % 1～2
                for j=1:Num_Tx % 1～2 Signal_disは後で改善する ⇒ ただデータを保存してるだけなので、Fadingに直接の影響はない。
                    [OutSignal1(:,:,j,u),Signal_dis1(j,u), Fading1(:,:,j,u)]=channel_only_AM(BaseSignal(:,:,j),NormMulti,MultiNo,...
                    Transrate,Doppler,Delay,SNR,com_delay,channel_rand(j,i,u,:,:,k,r));% channel_only([80,22,2],[7,1],7,1000,4,2,7,0)
                    OutSignal(:,:,i,u)=OutSignal(:,:,i,u)+OutSignal1(:,:,j,u); % [1,1760,2,4]
                    Signal_dis(i,u)=Signal_dis(i,u)+Signal_dis1(j,u); % [2,4]
                    Fading(:,:,i,u)=Fading(:,:,i,u)+Fading1(:,:,j,u); % [7,22,2,4]
                end
            end    
        end
        % OutSignal=[1,80*22,2]（マルチパス環境では、Multipathfading)
        % OutSignal[1,80*22,1]⇒h1*x1+h2'*x2
        % OutSignal[1,80*22,2]⇒h1'*x1+h2*x2
        % Signal_dis=送信信号[80,22,2]の標準偏差
        % Fading=[7,22,2] フェーシング　⇒ 各Pathのシンボルに対してのフェージング
        %%  AWGN channel
        TimeSignal=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        noise=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        for u=1:Num_User
            for receive=1:Num_Rx
                [TimeSignal(:,:,receive,u),noise(1,:,receive,u)]=awgn(OutSignal(:,:,receive,u),Signal_dis(receive,u),MultiNo,Doppler,SNR,Num_sym,...
                    NumCarr,guardtime,Num_pilot,coding_rate,perfect,wordsize);% [[1,22*80,2,4],1,7,4,7,20,64,16,2,cording_rate,0,1]
            end
        end

        % TimeSignal[1,80*22,1]⇒h1*x1+h2'*x2+noise1
        % TimeSignal[1,80*22,2]⇒h1'*x1+h2*x2+noise2
        clear OutSignal;% save memory
        clear BaseSignal;% save memory
  
         %% Rx device 
        [Datarx,Datarx_hd]=receiver_non_FSS_MU_MIMO_AM_main(TimeSignal,ifftsize,carriers,wordsize,guardtype,...
                                    guardtime,Num_sym,Num_pilot,Doppler,proc_gain,Num_Tx,Num_Rx,Num_User,noise,channel_ranking,user_index,block_num_per_user);% [[1,22*80,2,4],64,[1,64],1,2,4,20,2,4,64,2,2,[1,22*80,2,4]]
        
        % Datarx_hd=復調したデータ [20,64,2,2]
        % Datarx=実数部の生データを符号反転したもの [20,64,2,2]
        if strength_len==0
           Summary=error_count(Gen_data,Datarx_hd,Num_Tx);
           Result(k,1) = SNR;                  
            %% Add up the number of errorstea
           Result(k,l+1) = Result(k,l+1)+Summary(1);	
        else
           Summary=error_count(Gen_data,Datarx_hd,Num_Tx); % Gen_data=[20,64,2,2],Datarx_hd=[20,64,2,2],Num_Tx=2
           Summary2=error_count2(Ori_data,Datarx,g,int_pattern,Num_Tx); % Ori_data=[1,20*64*2-6],Datarx=[20,64,2,2],g=[2,7],int_pattern=[1,20*64*2*2]←index
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
legend('FSS-MIMO(MMSE_MLD)','FSS-MIMO(FEC)');
xlabel('Eb/No [dB]');
ylabel('BER');
axis([0,30,10^(-5),10^(0)]);
grid;
hold on;