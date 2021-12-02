function Summary=error_count2(Gen_data,Datarx,g,int_pattern,Num_Tx)
% Ori_data=[1,20*64*2-6],Datarx=[20,64,2,2],g=[2,7],int_pattern=[1,20*64*2*2]←index
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
   Data_int=reshape(Datarx,[NumSymb,NumCarr*wordsize,Num_Tx]);
end
Data=reshape(Data_int,1,NumSymb*NumCarr*wordsize*Num_Tx);

%****************
% Deinterleaving
%****************
for xx=1:size(Data,2)
   Data_fin(1,int_pattern(xx))=Data(1,xx);
end
clear Data;%save memory
mem=2^(const_length-1);
%Data_fin_2=(-1).^(Data_fin');
Data_fin_2=(Data_fin');
RxData=(Decoder(Data_fin_2, g,mem))';
clear Data_fin;%save memory
clear Data_fin_2;%save memory



%******************
%Calculate the BER
%******************
Errors1=find(Gen_data-RxData);
NumErr1=length(Errors1);
NumData1=size(RxData,1)*size(RxData,2);%
BER = NumErr1/NumData1;
if NumErr1~=0
  PER=1; 
else
  PER=0; 
end
Summary=[BER,PER];

%[Datatx(2,:);Datarx(2,:)]