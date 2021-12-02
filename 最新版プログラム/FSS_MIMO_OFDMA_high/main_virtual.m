function [channel_ranking,user_index]=main_virtual(SNR,Doppler,com_delay,ifftsize,guardtype,guardtime,...
    MultiNo,Transrate,Delay,ww,Num_sym,NumCarr,Num_pilot,strength_len,...
    carriers,proc_gain,WordSizes,NumSizes,k,Num_Tx,Num_Rx,Num_User,channel_rand,r)

rand('state',sum(100*clock));%random probability

%**********************
% FEC function command ⇒ 誤り訂正部
%**********************
if strength_len==0
    g=0;
    coding_rate=1;% 符号化率：有用な符号の割合（冗長成分を除く）
    vertual_filename= 'result_uncoded_10Hz.txt';%File to store the results in
else
    g=g12(strength_len);% generate matrix
    coding_rate=1/2;% R ⇒ 入力1に対して、出力が2つ得られる。
    vertual_filename= 'result_coded_10Hz.txt';%File to store the results in
    vertual_filename1= 'result_uncoded_10Hz.txt';%File to store the results in
    vertual_filename2= 'result_throughput_10Hz.txt';%File to store the results in
end

%*****************************
% Various modulation schemes.
%*****************************
for l = 1:NumSizes % 1
    wordsize = WordSizes(l);
    %************
    % Tx device
    %************
    if strength_len==0
        vertual_Gen_data=MIMO_data(Num_sym,wordsize,NumCarr,Num_Tx);% MIMO_data(20,2,64,2)⇒ 最終的に[20,64,2,2]の0or1が返ってくる
    else
        % Original Version 
        %[Gen_data,Ori_data,int_pattern]=CMIMO_data(Num_sym,wordsize,NumCarr,g,Num_Tx);% convolutional code(R=1/2,K=7)
        % Gen_data=[20,64,2,2]←0or1,Ori_data=[1,20*64*2-6]←0or1,int_pattern=[1,20*64*2*2]←index
        
        % Fixed Version → 全ストリームに一気にインタリーブを掛けている
        [vertual_Gen_data,vertual_Ori_data,vertual_int_pattern]=CMIMO_data(Num_sym,wordsize,NumCarr,g,Num_Tx);
        % CMIMO_data(20,2,64,[1,7],2);
        % convolutional code(R=1/2,K=7)
        % vertual_Gen_data=[20,64,2,2]←0or1,vertual_Ori_data=[1,20*64-6,2]←0or1,vertual_int_pattern=[20*64*2,2]←index
    end
    vertual_Datatx=transmitter_non_FSS_MIMO_OFDMA(vertual_Gen_data,Num_sym,wordsize,NumCarr,guardtype,guardtime,...
        Num_pilot,proc_gain,Num_Tx); % transmitter([20,64,2,2],20,2,64,2,16,2,64,2) ⇒ [22,64,2]のデータ（パイロット信号 含）
    
    % 拡散成分を3次元方向に足し合わせているので、いつもと異なった求め方になっていることに注意する。
    % **********************************************
    %  Make a timewave using IFFT and normalization
    % **********************************************
    
    vertual_BaseSignal=zeros(NumCarr,Num_sym+Num_pilot,Num_Tx); % [64,22,2]
    for i=1:Num_Tx % 1～8
        vertual_BaseSignal(:,:,i)=(ifft(vertual_Datatx(:,:,i)'))*sqrt(NumCarr);% ⇒ [64,22,i]（Datatxで正規化の為に、sqrt(64)で割ったので、ifft部では、sqrt(64)を掛けている）
    end
    % 複素共役転置を取った後にifftをしている。
    % ifft,fftはifft(x,[],1)で列方向のifft(fft)
    %           ifft(x,[],2)で行方向のifft(fft)
    % ※何も指定しない場合は、ifft(x,[],1)と同じ扱いになる。
    % [64,22,8]
    %*******************
    % Add a Guard Period
    %*******************
    if guardtype~=0% 等しくない場合は1を返し、等しい場合は0を返す。
        if guardtype == 1%guard period is just a zero transmission
            vertual_BaseSignal=[zeros(guardtime,Num_sym+Num_pilot,Num_Tx);vertual_BaseSignal];
        elseif guardtype == 2% ⇒ 64サブキャリアの後方成分を前方に持ってくる
            vertual_EndSignal = size(vertual_BaseSignal,1);% Find the number of columns in the BaseSignal ⇒ 64
            vertual_BaseSignal=[vertual_BaseSignal((vertual_EndSignal-guardtime+1):vertual_EndSignal,:,:); vertual_BaseSignal];% [80,22,2] ⇒ ガードインターバルをつけた
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
    % とりあえず、ユーザ1のみを考慮する⇒ユーザ数を拡大する
    vertual_OutSignal1=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Tx,Num_User); % [1,22*80,2,4]
    vertual_OutSignal=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
    vertual_Signal_dis1=zeros(Num_Tx,Num_User); % [2,4]
    vertual_Signal_dis=zeros(Num_Rx,Num_User);% [2,4]
    vertual_Fading1=zeros(MultiNo,Num_sym+Num_pilot,Num_Tx,Num_User); % [7,22,2,4]
    vertual_Fading=zeros(MultiNo,Num_sym+Num_pilot,Num_Rx,Num_User); % [7,22,2,4]
    for u=1:Num_User
        for i=1:Num_Rx % 1～2
            for j=1:Num_Tx % 1～2 Signal_disは後で改善する ⇒ ただデータを保存してるだけなので、Fadingに直接の影響はない。
                [vertual_OutSignal1(:,:,j,u),vertual_Signal_dis1(j,u), vertual_Fading1(:,:,j,u)]=channel_only_MIMO_OFDMA(vertual_BaseSignal(:,:,j),NormMulti,MultiNo,...
                    Transrate,Doppler,Delay,SNR,com_delay,channel_rand(j,i,u,:,:,k,r));
                % channel_only([80,22,2],[7,1],7,1000,4,2,7,0)
                % channel_rand=[2,2,4,2,7,7,rep],[Num_Tx,Num_Rx,Num_User,2,MultiNo,k,rep]
                % k：1～7,MultiNoの部分はchannel_only内で1～7を割り当てる,rep=繰り返し回数
                vertual_OutSignal(:,:,i,u)=vertual_OutSignal(:,:,i,u)+vertual_OutSignal1(:,:,j,u);
                vertual_Signal_dis(i,u)=vertual_Signal_dis(i,u)+vertual_Signal_dis1(j,u);
                vertual_Fading(:,:,i,u)=vertual_Fading(:,:,i,u)+vertual_Fading1(:,:,j,u);
            end
        end
    end
    % OutSignal=[1,80*22,2]（マルチパス環境では、Multipathfading)
    % OutSignal[1,80*22,2]⇒h1*x1+h2'*x2
    % OutSignal[1,80*22,2]⇒h1'*x1+h2*x2
    % ※hはチャネル応答を表す
    % Signal_dis=送信信号[80,8,2]の標準偏差 ⇒ 2つの値が得られる.
    % Fading=[7,22,2] フェーシング　⇒ 各Pathのシンボルに対してのフェージング
    vertual_TimeSignal=vertual_OutSignal; % [1,22*80,2,4]
    
    clear OutSignal;% save memory
    clear BaseSignal;% save memory
    
    %**************
    % Rx device
    %**************
    % ここで得たいのは、チャネル推定の値のみ
    H_m_resp=receiver_non_FSS_MIMO_OFDMA_vertual(vertual_TimeSignal,ifftsize,carriers,wordsize,guardtype,...
        guardtime,Num_sym,Num_pilot,Doppler,proc_gain,Num_Tx,Num_Rx,Num_User);% [[1,8*80,2],64,[1,64],1,2,16,20,8,10,64,8,2,4]
    % [22,64,2,4]
    
    
    H_m=H_m_resp; % [2,2,64,4]
    % [1,:,64,1] → ユーザ1のアンテナ1へ向かうときの伝播チャネル推定情報　
    % [2,:,64,1] → ユーザ1のアンテナ2へ向かうときの伝搬チャネル推定情報
    % ユーザ1にとってほしい情報は、[1,1,64,1]の部分の情報であり、[1,2,64,1]の部分は結局、MMSECで分離されるので無視

    %************************
    %   channel ranking
    %************************
    user=zeros(Num_User,Num_User,Num_Tx); % [4,4,2]
    proc_gain_block=NumCarr/proc_gain;
    for i=1:Num_Rx % 1～2
        for u=1:Num_User % 1～4
            for fscb=1:proc_gain_block % 1～4
               user(u,fscb,i)=abs(sum(H_m_resp(i,i,(fscb-1)*proc_gain+1:fscb*proc_gain,u),3)); 
            end
        end
    end      
    % 各サブキャリアブロック毎に平均値を算出する
    user_mean=zeros(Num_User,Num_Rx); % [1,4,2]
    
    for i=1:Num_Rx
        for u=1:Num_User
           user_mean(u,i)=abs(mean(user(u,:,i)));  % [4,2]
        end
    end
    
    % 平均値が小さいユーザからサブキャリアブロックを割り当てていく。
    
    M_abs=zeros(Num_User,Num_Tx,Num_Tx); % [4,2,2]
    I_abs=zeros(Num_User,Num_Tx,Num_Tx); % [4,2,2]
    user_index=zeros(Num_User,Num_Tx); % [4,2]
    channel_ranking=zeros(Num_User,Num_Tx); % [4,2]
    
    for i=1:Num_Tx % 1～2
        for u=1:Num_User % 1～4
            [M_abs(u,1,i),I_abs(u,1,i)]=min(user_mean(:,i));
            user_index(u,i)=I_abs(u,1,i);
            [M_abs(u,2,i),I_abs(u,2,i)]=max(user(I_abs(u,1,i),:,i));
            channel_ranking(u,i)=I_abs(u,2,i);
            % 評価し終わったら、ゼロ又は大きな数で埋める
            user(I_abs(u,1,i),:,i)=zeros(1,Num_User);
            user(:,I_abs(u,2,i),i)=zeros(Num_User,1);
            user_mean(I_abs(u,1,i),i)=100;
        end
    end
end
