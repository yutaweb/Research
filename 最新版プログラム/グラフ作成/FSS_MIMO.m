% SNR定義
SNRMin=0;
SNRInc=5;
SNRMax=30;

% FECを使っていない場合
coding_rate=1;
% Resultにはループ平均値が入っている → 間違っているので、修正したもので上書きする
R_FMZ=[]; % FSS-MIMO(ZF)
R_FMM=[]; % FSS-MIMO(MMSE)
R_MZ=[]; % MIMO(ZF)
R_MM=[]; % MIMO(MMSE)
R_MML=[0.170889388020833,0.066503645833333,0.014764713541667,0.001994140625,0.0002228515625,0.000021158854167,0.000003059895833]; % MIMO(MLD)
% 色は適当につけているので、あとで変更する
semilogy(SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),R_FMZ(1,:)','g*--',...
         SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),R_FMM(1,:)','bx--',...
         SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),R_MZ(1,:)','m>--',...
         SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),R_MM(1,:)','r>--')
%          SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),R_MM(1,:)','b>--'
%      );
legend('FSS-MIMO(ZF)',...
       'FSS-MIMO(MMSE)',...
       'MIMO(ZF)',...
       'MIMO(MMSE)');
%        'MIMO(MLD)');
xlabel('Eb/No [dB]');
ylabel('BER');
axis([0,30,10^(-5),10^(0)]);
grid;