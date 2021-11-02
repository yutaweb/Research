function [channel_ranking,user_index]=main_virtual(SNR,Doppler,com_delay,ifftsize,guardtype,guardtime,...
    MultiNo,Transrate,Delay,ww,Num_sym,NumCarr,Num_pilot,strength_len,...
    carriers,proc_gain,WordSizes,NumSizes,k,Num_Tx,Num_Rx,Num_User,channel_rand,r)

rand('state',sum(100*clock));%random probability

%**********************
% FEC function command 
%**********************
if strength_len==0
    g=0;
    coding_rate=1;
    vertual_filename= 'result_uncoded_10Hz.txt';%File to store the results in
else
    g=g12(strength_len);% generate matrix
    coding_rate=1/2;
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
        vertual_Gen_data=MIMO_data(Num_sym,wordsize,NumCarr,Num_Tx);% MIMO_data(20,2,64,2)
    else
        % Original Version 
        %[Gen_data,Ori_data,int_pattern]=CMIMO_data(Num_sym,wordsize,NumCarr,g,Num_Tx);% convolutional code(R=1/2,K=7)
        % Gen_data=[20,64,2,2] 0or1,Ori_data=[1,20*64*2-6] 0or1,int_pattern=[1,20*64*2*2] index
        
        [vertual_Gen_data,vertual_Ori_data,vertual_int_pattern]=CMIMO_data(Num_sym,wordsize,NumCarr,g,Num_Tx);
        % CMIMO_data(20,2,64,[1,7],2);
        % convolutional code(R=1/2,K=7)
        % vertual_Gen_data=[20,64,2,2] 0or1,vertual_Ori_data=[1,20*64-6,2] 0or1,vertual_int_pattern=[20*64*2,2] index
    end
    vertual_Datatx=transmitter_non_FSS_MIMO_OFDMA(vertual_Gen_data,Num_sym,wordsize,NumCarr,guardtype,guardtime,...
        Num_pilot,proc_gain,Num_Tx); % transmitter([20,64,2,2],20,2,64,2,16,2,64,2) [24,64,2]
    
    % **********************************************
    %  Make a timewave using IFFT and normalization
    % **********************************************
    
    vertual_BaseSignal=zeros(NumCarr,Num_sym+Num_pilot,Num_Tx); % [64,22,2]
    for i=1:Num_Tx % 1 - 4
        vertual_BaseSignal(:,:,i)=(ifft(vertual_Datatx(:,:,i)'))*sqrt(NumCarr);% [64,22,i]
    end
    % [64,22,4]
    %*******************
    % Add a Guard Period
    %*******************
    if guardtype~=0
        if guardtype == 1%guard period is just a zero transmission
            vertual_BaseSignal=[zeros(guardtime,Num_sym+Num_pilot,Num_Tx);vertual_BaseSignal];
        elseif guardtype == 2
            vertual_EndSignal = size(vertual_BaseSignal,1);% Find the number of columns in the BaseSignal 64
            vertual_BaseSignal=[vertual_BaseSignal((vertual_EndSignal-guardtime+1):vertual_EndSignal,:,:); vertual_BaseSignal];% [80,22,2]
        end
    end
    %*****************
    %  Channel block
    %*****************
    %%  Rayleigh fading channel
    Multi = zeros(MultiNo,1);% [7,1]
    for w=1:MultiNo% 1 - 7
        Multi(w)=10^(0-(((w-1)*ww)/10));% 10^0?ｽ?10^(-6/10)
    end
    Norm=sum(Multi); NormMulti=Multi./(Norm);% Normalization
    vertual_OutSignal1=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Tx,Num_User); % [1,22*80,2,4]
    vertual_OutSignal=zeros(1,(Num_sym+Num_pilot)*(NumCarr+guardtime),Num_Rx,Num_User); % [1,22*80,2,4]
    vertual_Signal_dis1=zeros(Num_Tx,Num_User); % [2,4]
    vertual_Signal_dis=zeros(Num_Rx,Num_User);% [2,4]
    vertual_Fading1=zeros(MultiNo,Num_sym+Num_pilot,Num_Tx,Num_User); % [7,22,2,4]
    vertual_Fading=zeros(MultiNo,Num_sym+Num_pilot,Num_Rx,Num_User); % [7,22,2,4]
    for u=1:Num_User
        for i=1:Num_Rx % 1 - 2
            for j=1:Num_Tx % 1 - 2 Signal_dis
                [vertual_OutSignal1(:,:,j,u),vertual_Signal_dis1(j,u), vertual_Fading1(:,:,j,u)]=channel_only_MIMO_OFDMA(vertual_BaseSignal(:,:,j),NormMulti,MultiNo,...
                    Transrate,Doppler,Delay,SNR,com_delay,channel_rand(j,i,u,:,:,k,r));
                % channel_only([80,22,2],[7,1],7,1000,4,2,7,0)
                % channel_rand=[2,2,4,2,7,7,rep],[Num_Tx,Num_Rx,Num_User,2,MultiNo,k,rep]
                vertual_OutSignal(:,:,i,u)=vertual_OutSignal(:,:,i,u)+vertual_OutSignal1(:,:,j,u);
                vertual_Signal_dis(i,u)=vertual_Signal_dis(i,u)+vertual_Signal_dis1(j,u);
                vertual_Fading(:,:,i,u)=vertual_Fading(:,:,i,u)+vertual_Fading1(:,:,j,u);
            end
        end
    end
    % OutSignal=[1,80*22,2]
    % OutSignal[1,80*22,2]h1*x1+h2'*x2
    % OutSignal[1,80*22,2]h1'*x1+h2*x2
    % Signal_dis=[80,8,2]
    % Fading=[7,22,2]
    vertual_TimeSignal=vertual_OutSignal; % [1,22*80,2,4]
    
    clear OutSignal;% save memory
    clear BaseSignal;% save memory
    
    %**************
    % Rx device
    %**************
    H_m_resp=receiver_non_FSS_MIMO_OFDMA_vertual(vertual_TimeSignal,ifftsize,carriers,wordsize,guardtype,...
        guardtime,Num_sym,Num_pilot,Doppler,proc_gain,Num_Tx,Num_Rx,Num_User);% [[1,8*80,2],64,[1,64],1,2,16,20,8,10,64,8,2,4]
    % [22,64,2,4]
    
    H_m=H_m_resp; % [2,2,64,4]
    % [1,:,64,1]
    % [2,:,64,1]

    %************************
    %   channel ranking
    %************************
    user=zeros(Num_User,Num_User,Num_Tx); % [4,4,2]
    proc_gain_block=NumCarr/proc_gain;

    % 各サブキャリアブロック毎に合計値を算出している（adaptiveに使用する）
%     for i=1:Num_Rx % 1 - 2
%         for u=1:Num_User % 1 - 4
%             for fscb=1:proc_gain_block % 1 - 4
%                user(u,fscb,i)=abs(sum(H_m_resp(i,i,(fscb-1)*proc_gain+1:fscb*proc_gain,u),3)); 
%             end
%         end
%     end
%     

    % 各サブキャリアブロック毎に平均値を算出している
%     for i=1:Num_Rx % 1 - 2
%         for u=1:Num_User % 1 - 4
%             for fscb=1:proc_gain_block % 1 - 4
%                user(u,fscb,i)=abs(mean(H_m_resp(i,i,(fscb-1)*proc_gain+1:fscb*proc_gain,u),3)); 
%             end
%         end
%     end      
    
    % 各サブキャリアブロック毎に２乗平均を算出している（一番良い）（best-adaptiveに使用する）
    for i=1:Num_Rx % 1 - 2
        for u=1:Num_User % 1 - 4
            for fscb=1:proc_gain_block % 1 - 4
               user(u,fscb,i)=abs(mean(abs(H_m_resp(i,i,(fscb-1)*proc_gain+1:fscb*proc_gain,u)),3)); 
            end
        end
    end      
    
    user_mean=zeros(Num_User,Num_Rx); % [1,4,2]
    
    % 各ユーザ毎に平均を取る（一番良い）(adaptive,best-adaptive)
    for i=1:Num_Rx
        for u=1:Num_User
           user_mean(u,i)=abs(mean(user(u,:,i)));  % [4,2]
        end
    end
    
    % 各ユーザ毎に合計を取る
%     for i=1:Num_Rx
%         for u=1:Num_User
%            user_mean(u,i)=abs(sum(user(u,:,i)));  % [4,2]
%         end
%     end

    % 各ユーザごとに中央値をとる
%     for i=1:Num_Rx
%         for u=1:Num_User
%            user_mean(u,i)=abs(median(user(u,:,i)));  % [4,2]
%         end
%     end
    
    % 各ユーザごとに最小値をとる
%     for i=1:Num_Rx
%         for u=1:Num_User
%            user_mean(u,i)=abs(min(user(u,:,i)));  % [4,2]
%         end
%     end
    
    % 各ユーザごとに最大値をとる
%     for i=1:Num_Rx
%         for u=1:Num_User
%            user_mean(u,i)=abs(max(user(u,:,i)));  % [4,2]
%         end
%     end
    
    M_abs=zeros(Num_User,Num_Tx,Num_Tx); % [4,2,2]
    I_abs=zeros(Num_User,Num_Tx,Num_Tx); % [4,2,2]
    user_index=zeros(Num_User,Num_Tx); % [4,2]
    channel_ranking=zeros(Num_User,Num_Tx); % [4,2]
    
    % 全ユーザにおいて平均値が小さいものからサブキャリアの割り当てを決定する
    for i=1:Num_Tx % 1 - 2
        for u=1:Num_User % 1 - 4
            [M_abs(u,1,i),I_abs(u,1,i)]=min(user_mean(:,i));
            user_index(u,i)=I_abs(u,1,i);
            [M_abs(u,2,i),I_abs(u,2,i)]=max(user(I_abs(u,1,i),:,i));
            channel_ranking(u,i)=I_abs(u,2,i);
            user(I_abs(u,1,i),:,i)=zeros(1,Num_User);
            user(:,I_abs(u,2,i),i)=zeros(Num_User,1);
            user_mean(I_abs(u,1,i),i)=100;
        end
    end
    
%     % 全ユーザにおいて平均値が大きいものからサブキャリアの割り当てを決定する
%     for i=1:Num_Tx % 1 - 2
%         for u=1:Num_User % 1 - 4
%             [M_abs(u,1,i),I_abs(u,1,i)]=max(user_mean(:,i));
%             user_index(u,i)=I_abs(u,1,i);
%             [M_abs(u,2,i),I_abs(u,2,i)]=max(user(I_abs(u,1,i),:,i));
%             channel_ranking(u,i)=I_abs(u,2,i);
%             user(I_abs(u,1,i),:,i)=zeros(1,Num_User);
%             user(:,I_abs(u,2,i),i)=zeros(Num_User,1);
%             user_mean(I_abs(u,1,i),i)=0;
%         end
%     end

    % user_meanでなく、user_sumで評価してみる
    % user_min,user_mediam,user_absも試してみる
    % 全ユーザの平均値の中央値からサブキャリアを割り当てる（最初は大きいほうへ）
    % 全ユーザの平均値の中央値からサブキャリアを割り当てる（最初は小さいほうへ）
end
