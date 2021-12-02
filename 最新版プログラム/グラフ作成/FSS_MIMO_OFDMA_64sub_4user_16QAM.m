% SNR定義
SNRMin=0;
SNRInc=5;
SNRMax=30;

% FECを使っていない場合
coding_rate=1;
% Resultにはループ平均値が入っている → 間違っているので、修正したもので上書きする
R_FMOPL=[0.2815990234375,0.17216748046875,0.07466533203125,0.0198046875,0.003396875,0.00059228515625,0.00013388671875]; % Proposal(Low)
R_FMOCL=[0.2667716796875,0.15997236328125,0.0756833984375,0.02901484375,0.009453515625,0.00313232421875,0.001030078125]; % Conventional(Low)
R_FMOCH=[0.225545703125,0.1244619140625,0.05521044921875,0.0197951171875,0.0063009765625,0.00207216796875,0.00073818359375]; % Conventional(High)
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