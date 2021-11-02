function [H_m_resp,beam_noise]=main_virtual_MU_MIMO(SNR,Doppler,com_delay,ifftsize,guardtype,guardtime,...
    MultiNo,Transrate,Delay,ww,Num_sym,NumCarr,Num_pilot,strength_len,...
    carriers,proc_gain,WordSizes,NumSizes,k,Num_Tx,Num_Rx,Num_User,channel_rand,r,noise_rand)

%*****************************
% Various modulation schemes.
%*****************************
for l = 1:NumSizes % 1
    wordsize = WordSizes(l);
    %************
    % Tx device
    %************
    if strength_len==0
        vertual_Gen_data=MIMO_data(Num_sym,wordsize,NumCarr,Num_Tx);% MIMO_data(20,2,64,4)⇒ 最終的に[20,64,2,4]の0or1が返ってくる
    else
        % ※要修正（テスト信号にはFECは不要）
        [vertual_Gen_data,vertual_Ori_data,vertual_int_pattern]=CMIMO_data(Num_sym,wordsize,NumCarr,g,Num_Tx,Num_User);
        % 各ユーザごとにインタリーブを掛けている
        % Gen_data[20,64,2,1:2] → ユーザ1に送信したい信号成分
        % Gen_data[20,64,2,3:4] → ユーザ2に送信したい信号成分
    end
    vertual_Datatx=transmitter_non_FSS(vertual_Gen_data,Num_sym,wordsize,NumCarr,guardtype,guardtime,...
        Num_pilot,proc_gain,Num_Tx,Num_User); 
    % transmitter([20,64,2,4],20,2,64,2,16,2,64,4) ⇒ [24,64,4]のデータ（パイロット信号 含）

    % 拡散成分を3次元方向に足し合わせているので、いつもと異なった求め方になっていることに注意する。
    % **********************************************
    %  Make a timewave using IFFT and normalization
    % **********************************************
    
    vertual_BaseSignal=zeros(NumCarr,Num_sym+Num_pilot,Num_Tx); % [64,24,4]
    for i=1:Num_Tx % 1～4
        vertual_BaseSignal(:,:,i)=(ifft(vertual_Datatx(:,:,i)'))*sqrt(NumCarr);% ⇒ [64,24,i]（Datatxで正規化の為に、sqrt(64)で割ったので、ifft部では、sqrt(64)を掛けている）
    end
    % 複素共役転置を取った後にifftをしている。
    % ifft,fftはifft(x,[],1)で列方向のifft(fft)
    %           ifft(x,[],2)で行方向のifft(fft)
    % ※何も指定しない場合は、ifft(x,[],1)と同じ扱いになる。
    %*******************
    % Add a Guard Period
    %*******************
    if guardtype~=0% 等しくない場合は1を返し、等しい場合は0を返す。
        if guardtype == 1%guard period is just a zero transmission
            vertual_BaseSignal=[zeros(guardtime,Num_sym+Num_pilot,Num_Tx);vertual_BaseSignal];
        elseif guardtype == 2% ⇒ 64サブキャリアの後方成分を前方に持ってくる
            vertual_EndSignal = size(vertual_BaseSignal,1);% Find the number of columns in the BaseSignal ⇒ 64
            vertual_BaseSignal=[vertual_BaseSignal((vertual_EndSignal-guardtime+1):vertual_EndSignal,:,:); vertual_BaseSignal];% [80,24,4] ⇒ ガードインターバルをつけた
        end
    end
    % vertual_BaseSignal[80,24,1:2] → ユーザ1
    % vertual_BaseSignal[80,24,3:4] → ユーザ2
    % Beamformingによってストリームを制限
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
    vertual_OutSignal1=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Tx,Num_User); % [1,24*80,4,2]
    vertual_OutSignal=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,24*80,2,2]
    vertual_Signal_dis1=zeros(Num_Tx,Num_User); % [4,2]
    vertual_Signal_dis=zeros(Num_Rx,Num_User);% [2,2]
    vertual_Fading1=zeros(MultiNo,Num_sym+Num_pilot,Num_Tx,Num_User); % [7,24,4,2]
    vertual_Fading=zeros(MultiNo,Num_sym+Num_pilot,Num_Rx,Num_User); % [7,24,2,2]
    for u=1:Num_User % 1～2
        for i=1:Num_Rx % 1～2
            for j=1:Num_Tx % 1～4 Signal_disは後で改善する ⇒ ただデータを保存してるだけなので、Fadingに直接の影響はない。
                [vertual_OutSignal1(:,:,j,u),vertual_Signal_dis1(j,u), vertual_Fading1(:,:,j,u)]=channel_only_AM(vertual_BaseSignal(:,:,j),NormMulti,MultiNo,...
                    Transrate,Doppler,Delay,SNR,com_delay,channel_rand(j,i,u,:,:,k,r));
                % channel_only([80,24,4],[7,1],7,1000,4,2,7,0)
                % channel_rand=[4,2,2,2,7,7,rep],[Num_Tx,Num_Rx,Num_User,2,MultiNo,k,rep]
                % k：1～7,MultiNoの部分はchannel_only内で1～7を割り当てる,rep=繰り返し回数
                vertual_OutSignal(:,:,i,u)=vertual_OutSignal(:,:,i,u)+vertual_OutSignal1(:,:,j,u);
                vertual_Signal_dis(i,u)=vertual_Signal_dis(i,u)+vertual_Signal_dis1(j,u);
                vertual_Fading(:,:,i,u)=vertual_Fading(:,:,i,u)+vertual_Fading1(:,:,j,u);
            end
        end
    end
    
    fake_TimeSignal=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,24*80,2,2]
    beam_noise=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,24*80,2,2]
    % 全ユーザの受信アンテナ数の合計は8本なので、それぞれに異なるノイズを加える必要がある
    for user=1:Num_User
        for receive=1:Num_Rx
            [fake_TimeSignal(:,:,receive,user),beam_noise(1,:,receive,user)]=awgn(vertual_OutSignal(:,:,receive,user),vertual_Signal_dis(receive,user),MultiNo,Doppler,SNR,Num_sym,...
                NumCarr,guardtime,Num_pilot,1,0,wordsize,k,r,user,receive,noise_rand);% [[1,22*80,2,4],1,7,4,7,20,64,16,2,cording_rate,0,1]
        end
    end
    % vertual_OutSignal=[1,80*24,2,2]（マルチパス環境では、Multipathfading)
    % ※hはチャネル応答を表す
    % vertual_Signal_dis=送信信号[80,8,2]の標準偏差 ⇒ 2つの値が得られる.
    % vertual_Fading=[7,24,2] フェーシング　⇒ 各Pathのシンボルに対してのフェージング
    vertual_TimeSignal=vertual_OutSignal; % [1,24*80,2,2]
    
    clear OutSignal;% save memory
    clear BaseSignal;% save memory
    
    %**************
    % Rx device
    %**************
    % ここで得たいのは、チャネル推定の値のみなので、ノイズは加えない
    % ユーザ間でチャネル状態の情報共有はできないので、基地局に戻した後に復元する
    H_m_resp=receiver_non_FSS_MU_MIMO_AM(vertual_TimeSignal,ifftsize,carriers,wordsize,guardtype,...
        guardtime,Num_sym,Num_pilot,Doppler,proc_gain,Num_Tx,Num_Rx,Num_User);% [[1,24*80,2,2],64,[1,64],1,2,16,20,8,10,64,8,2,4]
    % [4,4,64]
end
