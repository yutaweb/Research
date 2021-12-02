function Summary=error_count2(Gen_data,Datarx,g,int_pattern)
% Gen_data= [1,20*64-6]の0or1,Datarx=[20,64,2] ⇒ 0or1に復調していない実数成分（符号反転している）,g=[2,7],int_pattern=[1,20*64*2]のindex
% Gen_data:送信機でencoderを通してFEC処理をし、ソートしたデータを生成する前のデータ
% Datarx:受信機で全ての処理が終了したデータの実数成分（符号反転している）※0or1データに直していないデータ
% int_pattern:送信機でencoderを通してFEC処理をし、ソートするために[20,64,2]の0～1ランダムデータを作成してソートした際のindex
%**************************************************************
%Error_count calculates BER
%
% Copyright (c) Chang-Jun Ahn 2000 (jun@sasase.ics.keio.ac.jp)
%**************************************************************
wordsize=size(Datarx,3);% 2
NumCarr=size(Datarx,2);% 64
NumSymb=size(Datarx,1);% 20
const_length=size(g,2);% 7

%*****************************
%Decoding the received signal
%*****************************
if wordsize==1
   Data_int=Datarx;
elseif wordsize~=1
   Data_int=reshape(Datarx,[NumSymb,NumCarr*wordsize,1]);% [20,64*2]
end
Data=reshape(Data_int,1,NumSymb*NumCarr*wordsize);% [1,20*64*2] ⇒ 受信機ですべての処理が終了したデータの実数部反転データ

%****************
% Deinterleaving ⇒ インターリーブ：ランダムにデータを並べ替えて雑音への耐性を付与する。デインターリーブはその逆。
%****************
for xx=1:size(Data,2)% 1～20*64*2
   Data_fin(1,int_pattern(xx))=Data(1,xx); 
   % COFDM_dataで並べ替えたデータを元の位置(COFDM_dataにおけるData)に戻している。[1,20*64*2]
   % COFDM_dataにおけるデータはencoder.mを通しているので、FEC処理をする前のデータに戻すにはDecoderに通す必要がある。
end
clear Data;%save memory
mem=2^(const_length-1);% 2^6=64 
%Data_fin_2=(-1).^(Data_fin');
Data_fin_2=(Data_fin');%[20*64*2,1]※複素共役転置をしているが、Data_finは実数成分のみなので、転置をしていると考えればよい。
RxData=(Decoder(Data_fin_2, g,mem))';%Data_fin=[20*64*2,1],g=[2,7],mem=64 ⇒ 得られたデータを転置している。[20*64-6]
clear Data_fin;%save memory
clear Data_fin_2;%save memory



%******************
%Calculate the BER
%******************
Errors1=find(Gen_data-RxData);
NumErr1=length(Errors1);
NumData1=size(RxData,1)*size(RxData,2);%
BER = NumErr1/NumData1;
if NumErr1~=0 % エラーの数が0でない場合⇒エラーがある場合
  PER=1; 
else
  PER=0; 
end
Summary=[BER,PER];

%[Datatx(2,:);Datarx(2,:)]