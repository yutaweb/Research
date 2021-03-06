%****************************
%ã?ä½æ?è???¼é å±±è£å¤ª
%  ä¿®æ­£æ¥?¼?2020/12/02
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
MultiNo=7;% ãã«ããã¹æ°
Transrate=20*10^6;
Delay=2;%Number of delay sample
ww=1;%dB degradation â? dBå£å?

NumCarr=64;
Num_sym=20;%Number of data symbols    
Num_pilot=2;%Number of pilot symbols
   
%% FEC Function Command

strength_len=0;% ææé·
if strength_len==0
   g=0;
   coding_rate=1;% ç¬¦å·åç?¼æç¨ãªç¬¦å·ã®å²åï¼å?é·æå?ãé¤ãï¼?
   filename= 'result_uncoded_10Hz.txt';%File to store the results in
else
   g=g12(strength_len);% generate matrix
   coding_rate=1/2;% R â? å¥å?1ã«å¯¾ãã¦ãå?ºåã2ã¤å¾ãããã?
   filename= 'result_coded_10Hz.txt';%File to store the results in
   filename1= 'result_uncoded_10Hz.txt';%File to store the results in
   filename2= 'result_throughput_10Hz.txt';%File to store the results in
end

%% Carriers used for a single wide OFDM channel

StartCarr=1;
FinCarr=ifftsize;% FinCarr=64
carriers=[StartCarr:CarrSpacing:FinCarr];% 0?½?64ã¾ã§ã?1ãééã§carriersã«ã?ã¼ã¿ãå?¥ã? [1,64]
rep=1000;%Number of loop times âã??ä»»æã?®åæ°ã«ã¼ãããï¼ã·ãã¥ã¬ã¼ã·ã§ã³åæ°?¼?

SNRMin=0;
SNRMax=30;
SNRInc=5;
perfect=0;%perfect channel option(1: perfect,0:non perfect)
proc_gain=8; %OFDMAã«ãã£ã¦ã?4ã¦ã¼ã¶ã«å²ãå½ã¦ãã?®ã§16Ã?16ã®ãããªã?ã¯ã¹
WordSizes=[2];%%BPSK?¼?1,QPSK?¼?2,16QAM?¼?4
NumSizes=length(WordSizes);% Numsizes=1 âã??è¡å?ã?®å¤§ãã
Num_Tx=2;%Number of transmiter â? éä¿¡ã¢ã³ã?ãã?®æ° 
Num_Rx=2;%Number of receiver â? åä¿¡ã¢ã³ã?ãã?®æ°
Num_User=8;%Number of user â? ã¦ã¼ã¶ã®æ°

% æ¡ä»¶ãçµ±ä¸?ãããããã«åã«ãã¹ã¦ã®ã©ã³ã?ã?å¤ãçæãã?
rng('default');
rng(10); % ã©ã³ã?ã?ã®å¤ãå¤æ´ããã?æã?¯ããã®æ°å­ãå¤æ´ãã
channel_rand=rand([Num_Tx,Num_Rx,Num_User,2,MultiNo,(SNRMax-SNRMin)/SNRInc+1,rep]); 
% [2,2,4,7,7,10] â? [éä¿¡ã¢ã³ã?ãæ°,åä¿¡ã¢ã³ã?ãæ°,ã¦ã¼ã¶ã¼æ°,channel_onlyå?ã§ã®ä½¿ç¨æ°,ãã£ãã«ä¼æ­è·¯æ°,ãã«ããã¹æ°,SNRåå¾©æ°,å¨åå¾©æ°]

headerstr='SNR (samples)';      
BERstr=[];
PhErrstr=[];

Result = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)
Result1 = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)
Result2 = zeros(floor((SNRMax-SNRMin)/SNRInc)+1,NumSizes+1);% zeros(7,2)

%% Various modulation schemes.

for l = 1:NumSizes % 1
   wordsize = WordSizes(l);% wordsize = 2
   disp(['wordsize: ' num2str(wordsize)]);% wordsize : 2ãè¡¨ç¤ºããã?
   BERstr = [BERstr ', BER b/Hz: ', int2str(wordsize)];% int2strã¯æ´æ°ãæå­ã«å¤æãã
   PhErrstr = [PhErrstr ', Ph b/Hz: ', int2str(wordsize)];   
  
   %% Repeat various SNR situation
  
   for k = 1:((SNRMax-SNRMin)/SNRInc)+1 % 1?½?7 
      SNR = (k-1)*SNRInc + SNRMin;% 0?½?30 â»SNRInc=3
	  disp(['   Eb/N0: ' num2str(SNR)]);% 0?½?30ã¾ã§ãè¡¨ç¤ºãã
      refData=[]; symwaves=[]; 
      
      %% Repeat each run several times
      for r = 1:rep % 1?½?100 â? å?ã?ã®Eb/Noã«å¯¾ãã¦ä»»æã?®åæ°ã·ãã¥ã¬ã¼ã·ã§ã³ãã
        disp(['      Repeat ', num2str(r)]);
        %% Tx device
        if strength_len==0
           Gen_data=MIMO_data(Num_sym,wordsize,NumCarr,Num_Tx);% MIMO_data(20,2,64,2)â? æ?çµçã«[20,64,2,2]ã®0or1ãè¿ã£ã¦ãã
        else
           [Gen_data,Ori_data,int_pattern]=CMIMO_data(Num_sym,wordsize,NumCarr,g,Num_Tx);% convolutional code(R=1/2,K=7)
           % Gen_data=[20,64,2,2]â?0or1,Ori_data=[1,20*64*2-6]â?0or1,int_pattern=[1,20*64*2*2]âindex
           % FECåã?³ã¤ã³ã¿ãªã¼ããæããéã¯ãMU_MIMOã®å ´åã?å?¨ã¹ããªã¼ã?ã«åããã?®ãæãã¦ãããã?®ã¨æãããã?
        end
        
        %% estimate channel
        [channel_ranking,user_index]=main_virtual(SNR,Doppler,com_delay,ifftsize,guardtype,guardtime,...
            MultiNo,Transrate,Delay,ww,Num_sym,NumCarr,Num_pilot,strength_len,carriers,proc_gain,...
            WordSizes,NumSizes,k,Num_Tx,Num_Rx,Num_User,channel_rand,r);
         % best_adaptiveðßé½ßÉp^[IÉ»è·éB
        %% transmitter
         Datatx=transmitter_non_FSS_MU_MIMO(Gen_data,Num_sym,wordsize,NumCarr,guardtype,guardtime,...
                          Num_pilot,proc_gain,Num_Tx); % transmitter([20,64,2,2],20,2,64,2,16,2,64,4) â? [22,64,2]ã®ã?ã¼ã¿?¼ãã¤ã­ã?ãä¿¡å· å«?¼?
        
        % å?ã¹ããªã¼ã?ã«ããã¦64ãµãã­ã£ãªã¢ãã¦ã¼ã¶æ°å?ã?ãå?å²ãã¦ãåã¦ã¼ã¶ã®ã?ã¼ã¿ã«å¯¾ãã¦æ¡æ£ç¬¦å·ãå²ãå½ã¦ã?
        % [22,64,2]
        % 4ã¦ã¼ã¶ã®å ´å?
        % [22,1 ?½?16,2] â? ã¦ã¼ã¶1ã«å²ãå½ã¦ã?
        % [22,17?½?32,2] â? ã¦ã¼ã¶2ã«å²ãå½ã¦ã?
        % [22,33?½?48,2] â? ã¦ã¼ã¶3ã«å²ãå½ã¦ã?
        % [22,49?½?64,2] â? ã¦ã¼ã¶4ã«å²ãå½ã¦ã?
        % 8ã¦ã¼ã¶ã®å ´åã?¯ãµãã­ã£ãªã¢ã?8å?å²ã?16ã¦ã¼ã¶ã®å ´åã?¯ãµãã­ã£ãªã¢ã?16å?å²ããã?
        % â»ããã§ãæ³¨æããªããã°ãªããªã?ã®ã¯ãIEEEã®è¦æ?¼ã§ã¯ããµãã­ã£ãªã¢éã«ãã«ãè¨­ããå¿?è¦ãããã®ã§ãå®éã«ä½¿ããã®ã¯64ãµãã­ã£ãªã¢ã®å ´åã?¯ã?52å?
        % æ¡æ£æå?ã?3æ¬¡å?æ¹åã«è¶³ãåããã¦ã?ãã?®ã§ãå¾æ¥ã¨ç°ãªã£ãæ±ãæ¹ã«ãªã£ã¦ã?ããã¨ã«æ³¨æããã??
        %%  Make a timewave using IFFT and normalization
        % å¨ãµãã­ã£ãªã¢ã«ä¸?æ°ã«IFFTãæããå ´å? M
        BaseSignal=zeros(NumCarr,Num_sym+Num_pilot,Num_Tx); % [64,22,2]
        for i=1:Num_Tx % 1?½?2
            BaseSignal(:,:,i)=(ifft(Datatx(:,:,i)'))*sqrt(NumCarr);% â? [64,22,i]?¼?Datatxã§æ­£è¦åã®çºã«ãsqrt(64)ã§å²ã£ãã?®ã§ãiffté¨ã§ã¯ãsqrt(64)ãæãã¦ã?ãï¼?
        end
         
        % å?ã¦ã¼ã¶ã«å²ãå½ã¦ããµãã­ã£ãªã¢ãã­ã?ã¯ãã¨ã«IFFTãæããå ´å?
%         proc_gain_block=NumCarr/proc_gain; % 4
%         for fscb=1:proc_gain_block% 1?½?4
%             for i=1:Num_Tx % 1?½?2
%                 BaseSignal((fscb-1)*proc_gain+1:fscb*proc_gain,:,i)=(ifft(Datatx(:,(fscb-1)*proc_gain+1:fscb*proc_gain,i)'))*sqrt(proc_gain);
%             end
%         end
     
        % è¤?ç´?å±å½¹è»¢ç½®ãåã£ãå¾ã«ifftããã¦ã?ãã??
        % ifft,fftã¯ifft(x,[],1)ã§åæ¹åã?®ifft(fft)
        %           ifft(x,[],2)ã§è¡æ¹åã?®ifft(fft)
        % â»ä½ãæ?å®ããªã?å ´åã?¯ãifft(x,[],1)ã¨åãæ±ã?ã«ãªãã??
        % [64,22,2]
        
        %% Add a Guard Period
		if guardtype~=0% ç­ãããªã?å ´åã?¯1ãè¿ããç­ãã?å ´åã?¯0ãè¿ãã?
            if guardtype == 1%guard period is just a zero transmission
              BaseSignal=[zeros(guardtime,Num_sym+Num_pilot,Num_Tx);BaseSignal];
            elseif guardtype == 2% â? 64ãµãã­ã£ãªã¢ã®å¾æ¹æå?ãåæ¹ã«æã£ã¦ãã
      		  EndSignal = size(BaseSignal,1);% Find the number of columns in the BaseSignal â? 64
      		  BaseSignal=[BaseSignal((EndSignal-guardtime+1):EndSignal,:,:); BaseSignal];% [80,22,2] â? ã¬ã¼ãã¤ã³ã¿ã¼ãã«ãã¤ãã
            end	
        end        
        %%  Rayleigh fading channelâãã©ã?ããã§ã¼ã¸ã³ã°?¼æéã«ã¯ä¾å­ããªã?&å¨æ³¢æ°ã«å¯¾ãã¦ä¸?æ§ï¼?
        Multi = zeros(MultiNo,1);% [7,1]ã®0
        for w=1:MultiNo% 1?½?7
           Multi(w)=10^(0-(((w-1)*ww)/10));% 10^0?½?10^(-6/10)âdBãç´ãã
        end
        Norm=sum(Multi); NormMulti=Multi./(Norm);% Normalizationâæ­£è¦åãã¦ã?ãï¼?Multiã®åè¨å?¤ã?1ã«ãªãï¼?
        % ä»¥ä¸ã?ã¦ã¼ã¶æ°å?ã?ãä½ã
        OutSignal1=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        OutSignal=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        Signal_dis1=zeros(Num_Rx,Num_User); % [2,4]
        Signal_dis=zeros(Num_Rx,Num_User); % [2,4]
        Fading1=zeros(MultiNo,Num_sym+Num_pilot,Num_Rx,Num_User); % [7,22,2,4]
        Fading=zeros(MultiNo,Num_sym+Num_pilot,Num_Rx,Num_User); % [7,22,2,4]
        for u=1:Num_User % 1?½?4 â? å?ã¦ã¼ã¶ãã¨ã«ç°ãªãä¼æ¬ç°å¢?åã?³åä¿¡ã¢ã³ã?ããã¤ãºãæã£ã¦ã?ããã¨ã«èµ·å?ãã
            for i=1:Num_Rx % 1?½?2
                for j=1:Num_Tx % 1?½?2 Signal_disã¯å¾ã§æ¹å?ãã â? ãã ã?ã¼ã¿ãä¿å­ãã¦ãã ããªã®ã§ãFadingã«ç´æ¥ã®å½±é¿ã¯ãªã?ã?
                    [OutSignal1(:,:,j,u),Signal_dis1(j,u), Fading1(:,:,j,u)]=channel_only_AM(BaseSignal(:,:,j),NormMulti,MultiNo,...
                    Transrate,Doppler,Delay,SNR,com_delay,channel_rand(j,i,u,:,:,k,r));% channel_only([80,22,2],[7,1],7,1000,4,2,7,0)
                    OutSignal(:,:,i,u)=OutSignal(:,:,i,u)+OutSignal1(:,:,j,u); % [1,1760,2,4]
                    Signal_dis(i,u)=Signal_dis(i,u)+Signal_dis1(j,u); % [2,4]
                    Fading(:,:,i,u)=Fading(:,:,i,u)+Fading1(:,:,j,u); % [7,22,2,4]
                end
            end    
        end
        % OutSignal=[1,80*22,2]?¼ã?ã«ããã¹ç°å¢?ã§ã¯ãMultipathfading)
        % OutSignal[1,80*22,1]âh1*x1+h2'*x2
        % OutSignal[1,80*22,2]âh1'*x1+h2*x2
        % Signal_dis=éä¿¡ä¿¡å·[80,22,2]ã®æ¨æºåå·®
        % Fading=[7,22,2] ãã§ã¼ã·ã³ã°ã?â? åPathã®ã·ã³ãã«ã«å¯¾ãã¦ã®ãã§ã¼ã¸ã³ã°
        %%  AWGN channel
        TimeSignal=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        noise=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
        for u=1:Num_User
            for receive=1:Num_Rx
                [TimeSignal(:,:,receive,u),noise(1,:,receive,u)]=awgn(OutSignal(:,:,receive,u),Signal_dis(receive,u),MultiNo,Doppler,SNR,Num_sym,...
                    NumCarr,guardtime,Num_pilot,coding_rate,perfect,wordsize);% [[1,22*80,2,4],1,7,4,7,20,64,16,2,cording_rate,0,1]
            end
        end

        % TimeSignal[1,80*22,1]âh1*x1+h2'*x2+noise1
        % TimeSignal[1,80*22,2]âh1'*x1+h2*x2+noise2
        clear OutSignal;% save memory
        clear BaseSignal;% save memory
  
         %% Rx device 
        [Datarx,Datarx_hd]=receiver_non_FSS_MU_MIMO_AM_main(TimeSignal,ifftsize,carriers,wordsize,guardtype,...
                                    guardtime,Num_sym,Num_pilot,Doppler,proc_gain,Num_Tx,Num_Rx,Num_User,noise,channel_ranking,user_index);% [[1,22*80,2,4],64,[1,64],1,2,4,20,2,4,64,2,2,[1,22*80,2,4]]
        
        % Datarx_hd=å¾©èª¿ããã?ã¼ã¿ [20,64,2,2]
        % Datarx=å®æ°é¨ã®çãã¼ã¿ãç¬¦å·åè»¢ãããã?® [20,64,2,2]
        if strength_len==0
           Summary=error_count(Gen_data,Datarx_hd,Num_Tx);
           Result(k,1) = SNR;                  
            %% Add up the number of errors
           Result(k,l+1) = Result(k,l+1)+Summary(1);	
        else
           Summary=error_count(Gen_data,Datarx_hd,Num_Tx); % Gen_data=[20,64,2,2],Datarx_hd=[20,64,2,2],Num_Tx=2
           Summary2=error_count2(Ori_data,Datarx,g,int_pattern,Num_Tx); % Ori_data=[1,20*64*2-6],Datarx=[20,64,2,2],g=[2,7],int_pattern=[1,20*64*2*2]âindex
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