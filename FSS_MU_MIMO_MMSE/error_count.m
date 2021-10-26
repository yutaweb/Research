function BER=error_count(Gen_data,Datarx,Num_Tx)
% Gen_data=[20,64,2,2],Datarx=[20,64,2,2]
%**************************************************************
% Error_count calculates BER
%
% Copyright (c) Chang-Jun Ahn 2000 (jun@sasase.ics.keio.ac.jp)
%**************************************************************
wordsize=size(Datarx,3); % 2

%******************
%Calculate the BER
%******************
Errors1=find(Gen_data-Datarx); % 元のデータと比べた時の誤り
NumErr1=length(Errors1); % エラーの数
NumData1=size(Datarx,1)*size(Datarx,2)*wordsize*Num_Tx; % 20*64*2*2
BER = NumErr1/NumData1; % BERの値


%[Gen_data(1,:,:); Datarx(1,:,:)]

