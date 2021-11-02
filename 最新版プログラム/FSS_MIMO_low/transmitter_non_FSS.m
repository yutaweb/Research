function Datatx=transmitter_non_FSS(Gen_data,Numsymb,wordsize,NumCarr,guardtype,...
                                guardtime,Num_pilot,proc_gain,Num_Tx)
% Gen_data=[20,64,2,2]の0or1,Numsymb=20,wordsize=2,NumCarr=64,guardtype=2,guardtime=16,Num_pilot=2,proc_gain=64,Num_Tx=2
% Gen_data=[20,64,2,1]はアンテナ1からの成分
% Gen_data=[20,64,2,2]はアンテナ2からの成分
%*******************************************************************************         
% This function generates an OFDM time waveform based on the input parameters
% and the data given. The data is transmformed into a single frame OFDM signal. 
%
% INPUTS:
% ========
% Gen_data=Input data.
% Numsymb=Number of symbols.
% wordsize=modulation command. 
% NumCarr=Number of carriers.⇒搬送波の数
% guardtype=What type of guard period to use options:
%		      0: No Guard period
%		      1: zero level guard period
%		      2: cyclic extension of end of symbols ⇒ 信号の終端の巡回拡大
%              ⇒ 伝送波（サブキャリア）の最後の16個（仮）をコピーしてガードインターバルとして使う
% guardtime=Number of sample to use for the total guard time.
% Num_pilot=Number of pilot signals. ⇒ シンボルの前につける（シンボルサイズは20）
% Num_Tx=transmit anttena			
%
% OUTPUTS:
% ========
% Datatx=Output time signal for the OFDM waveform. 
%	  
% Copyright (c) 2000-2003 Chang-Jun Ahn (junny700@crl.go.jp)
%*******************************************************************************         
%****************
% Initialization
%****************
rand('state',sum(100*clock));% to increase the probability of random
proc_gain_block=NumCarr/proc_gain;% 64/64=1
spread=hadamard(proc_gain);% [64,64]のアダマール行列

%*********************
% Make a pilot signal
%*********************
Pilot=zeros(Num_pilot,NumCarr,Num_Tx);% [2,64,2]の1 ⇒ パイロット信号（受信機での信号特定に使用）
H=hadamard(Num_Tx);
for i=1:Num_Tx % 1～2
    Pilot(:,:,i)=repmat(H(:,i),1,NumCarr);
end
%hadamard(2)=[[1,1],[1,-1]]
%Pilot=[2,64,2]
%Pilot(:,:,1)=[1;1]が64個
%Pilot(:,:,2)=[1;-1]が64個

%***************************
% Modulation for input data
%***************************
Data=Gen_data*2-1;% [20,64,2,2]の-1or1
clear Gen_data;% Gen_dataの中身を空にする

%*********************
% Packet data mapping
%*********************
Datatx=zeros(Numsymb+Num_pilot,NumCarr,Num_Tx);% [22,64,2]
Datatx1=zeros(Numsymb,NumCarr,Num_Tx);% [20,64,2]
Datatx(1:Num_pilot,:,:)=Pilot;% [2,64,2]にパイロット信号を入れる
if wordsize==1% BPSK ⇒ [20,64,2] 2by2
   Datatx1=Data;% [20,64,2]の-1or1
elseif wordsize==2% QPSK ⇒ [20,64,2] 2by2 ※ Data = [20,64,2,2]
    for transmit=1:Num_Tx
        Datatx1(:,:,transmit)=(Data(:,:,1,transmit)+sqrt(-1).*Data(:,:,2,transmit))*sqrt(2); % [20,64,2]の-1-i,-1+i,1-i,1+i
    end
elseif wordsize==4% 16QAM
    for transmit=1:Num_Tx
        Datatx1(:,:,transmit)=((Data(:,:,1,transmit)*2+Data(:,:,2,transmit))+sqrt(-1).*(Data(:,:,3,transmit)*2+Data(:,:,4,transmit)))./sqrt(5);
    end
end
% Dataは[NumSymb,NumCarr,wordsize,Num_Tx]となっているので、wordsizeに該当する箇所を数字で置きなおす

% [20,64,64,2]を3次元方向に和を取る。⇒ [20,64,2] これを正規化するために、sqrt(64)で割る。
% [20,64,1]には送信アンテナ１からのシンボル、[20,64,2]には送信アンテナ2からのシンボル
Datatx(Num_pilot+1:Numsymb+Num_pilot,:,:)=Datatx1/Num_Tx;% [20,64,2]のデータ
% Datatx(22,64,2)のデータが返る。

% 周波数拡散したデータが入っているので、訳が分からなくなっている。

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[temp_1, int_pattern_1]=sort(rand(1,size(Datatx,2)));% [] sort：小さい順に並べる
%for xx=1:size(Datatx,2)
%   Datatxtx(:,xx)=Datatx(:,int_pattern_1(xx));
%end
% temp_1:値,int_pattern_1:インデックス







