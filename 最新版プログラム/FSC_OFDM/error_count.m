function BER=error_count(Gen_data,Datarx)
% Gen_data=[20,64,2],Datarx=[20,64,2]
%**************************************************************
%Error_count calculates BER
%
% Copyright (c) Chang-Jun Ahn 2000 (jun@sasase.ics.keio.ac.jp)
%**************************************************************
wordsize=size(Datarx,3);% 2

%******************
%Calculate the BER
%******************
Errors1=find(Gen_data-Datarx);% findは非ゼロの要素を取り出してindexを配列として格納する。
NumErr1=length(Errors1);% エラーの数を数える。
NumData1=size(Datarx,1)*size(Datarx,2)*wordsize;%
BER = NumErr1/NumData1;


%[Gen_data(1,:,:); Datarx(1,:,:)]

