% SNR定義
SNRMin=0;
SNRInc=5;
SNRMax=30;

% FECを使っていない場合
coding_rate=1;
% Resultにはループ平均値が入っている → 間違っているので、修正したもので上書きする
R_FMOPL=[]; % Proposal(Low)
R_FMOCL=[]; % Conventional(Low)
R_FMOCH=[]; % Conventional(High)
% 色は適当につけているので、あとで変更する
semilogy(SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),R_FMOPL(1,:)','g*--',...
         SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),R_FMOCL(1,:)','bx--',...
         SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),R_FMOCH(1,:)','m>--');
legend('Proposal(Low-granularity)',...
       'Conventional(Low-granularity)',...
       'Conventional(High-granularity)');
xlabel('Eb/No [dB]');
ylabel('BER');
axis([0,30,10^(-5),10^(0)]);
grid;