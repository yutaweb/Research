function Datatx=transmitter(Gen_data,Numsymb,wordsize,NumCarr,guardtype,...
                                guardtime,Num_pilot,proc_gain,Num_Tx)
% Gen_data=[20,64,2,2]の0or1,Numsymb=20,wordsize=2,NumCarr=64,guardtype=2,guardtime=16,Num_pilot=2,proc_gain=16,Num_Tx=2
% ユーザ数が4の場合は、サブキャリアを4分割して各ユーザに割り当てる
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
proc_gain_block=NumCarr/proc_gain;% 64/16=4
spread=hadamard(proc_gain);% [16,16]のアダマール行列 → 各ユーザに対して拡散符号をかける

%*********************
% Make a pilot signal
%*********************
Pilot=zeros(Num_pilot,NumCarr,Num_Tx);% [2,64,2]の0 ⇒ パイロット信号（受信機での信号特定に使用）
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
        Datatx1(:,:,transmit)=(Data(:,:,1,transmit)+sqrt(-1).*Data(:,:,2,transmit)); % [20,64,2]の-1-i,-1+i,1-i,1+i
    end
elseif wordsize==4% 16QAM
    for transmit=1:Num_Tx
        Datatx1(:,:,transmit)=((Data(:,:,1,transmit)*2+Data(:,:,2,transmit))+sqrt(-1).*(Data(:,:,3,transmit)*2+Data(:,:,4,transmit)))./sqrt(5);
    end
end

for transmit=1:Num_Tx
    for fscb=1:proc_gain_block% 1～4 ※ユーザ数が4
       for fsc=1:proc_gain% 1～16 ※ユーザ数が4
          Datatx2(:,(fscb-1)*proc_gain+1:fscb*proc_gain,fsc,transmit)=...
          repmat(Datatx1(:,(fscb-1)*proc_gain+fsc,transmit),1,proc_gain).*repmat(spread(fsc,:),Numsymb,1);
          % Datatx1には[20,64,1]の-1or1が入っているので、Datatx1(20,1～64)を[1,64]の配列になるようにコピーする。
          % 各サブキャリアの成分を64サブキャリア分複製して周波数拡散している
          % 64個の[20,64]⇒[20,64,1]の成分はDatatx1の[20,1]を64個複製したものとアダマール行列の1列目[1,64]を20個複製したものの乗算
          % spread(1～64,64)を[20,1]に配列なるようにコピーする。64個の[20,64]
          % OFDM：狭帯域を作り出すことでISIの影響を低減する。この時、伝送速度が低下するので、並列伝送によって解決する。
       end
    end
end
% 上記のループ処理について ※ユーザ数が4
% 夫々に4分割された拡散符号が掛けられている。
% fscb=1の時、[20,1 ～16,1～16,2]
% fscb=2の時、[20,17～32,1～16,2]
% fscb=3の時、[20,33～48,1～16,2]
% fscb=4の時、[20,49～64,1～16,2]

Datatx2=sum(Datatx2,3)./sqrt(proc_gain);% [20,64,1,2]
% 4ユーザの場合は以下のように割り当てる
% [20,1 ～16,1,2] ⇒ ユーザ1に割り当てる
% [20,17～32,1,2] ⇒ ユーザ2に割り当てる
% [20,33～48,1,2] ⇒ ユーザ3に割り当てる
% [20,49～64,1,2] ⇒ ユーザ4に割り当てる
% 8ユーザの場合は8分割、16ユーザの場合は16分割
% [20,64,64,2]を3次元方向に和を取る。⇒ [20,64,2] これを正規化するために、sqrt(64)で割る。
Datatx(Num_pilot+1:Numsymb+Num_pilot,:,:)=permute(Datatx2,[1 2 4 3]);% [22,64,2]のデータ
% Datatx(22,64,2)のデータが返る。
% 周波数拡散したデータが入っているので、訳が分からなくなっている。

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[temp_1, int_pattern_1]=sort(rand(1,size(Datatx,2)));% [] sort：小さい順に並べる
%for xx=1:size(Datatx,2)
%   Datatxtx(:,xx)=Datatx(:,int_pattern_1(xx));
%end
% temp_1:値,int_pattern_1:インデックス







