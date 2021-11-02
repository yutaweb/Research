%=========================================
% シミュレーション定義
% 対象：FSS-MU-MIMO
% リピート回数：5000
% 日付：12/31
% シミュレーター名：Development-sim03-2-high
%=========================================
% SNR定義
SNRMin=0;
SNRInc=5;
SNRMax=30;
% FECを使っていない場合
coding_rate=1;
% Resultにはループ平均値が入っている
R_MO=[0.135917128906251,0.06607419921875,0.02614501953125,0.0089596484375,0.0028872265625,0.00093087890625,0.0003016796875]; % MIMO-OFDMA(MMSEC)
R_FMON=[0.16315587890625,0.0722305078125,0.01913818359375,0.00290802734375,0.0003573046875,0.0000537109375,0.00001619140625]; % FSS-MIMO-OFDMA(non-adaptive)
R_FMOA=[0.152721328125,0.063981875,0.01557650390625,0.0021634375,0.0002616015625,0.00004578125,0.0000121875]; % FSS-MIMO-OFDMA(adaptive-low)
R_FMOAL=[0.15818482421875,0.0681793359375,0.01746041015625,0.0025930859375,0.00031255859375,0.00005142578125,0.00001064453125]; % FSS-MIMO-OFDMA(best-adaptive-low)
% 色は適当につけているので、あとで変更する
semilogy(SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),R_MO(1,:)','g*--',...
         SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),R_FMON(1,:)','bx--',...
         SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),R_FMOA(1,:)','m>--',...
         SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),R_FMOAL(1,:)','r>--');
legend('MIMO-OFDMA',...
       'FSS-MIMO-OFDMA(non-adaptive)',...
       'FSS-MIMO-OFDMA(adaptive-low)',...
       'FSS-MIMO-OFDMA(best-adaptive-low)');
xlabel('Eb/No [dB]');
ylabel('BER');
axis([0,30,10^(-5),10^(0)]);
grid;