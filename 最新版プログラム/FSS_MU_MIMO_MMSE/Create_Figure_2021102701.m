%=========================================
% シミュレーション定義
% 対象：FSS-MU-MIMO
% リピート回数：500
% 日付：10/27
%=========================================
% SNR定義
SNRMin=0;
SNRInc=5;
SNRMax=30;
% FECを使っていない場合
coding_rate=1;
% Resultにはループ平均値が入っている
R_MUMIMO=[0.439120117187501,0.3675572265625,0.2487953125,0.099353515625,0.0100740234375,0.0003091796875,0.000075];
R_FSSMUMIMO=[0.4391197265625,0.367764453125,0.2493248046875,0.0991751953125,0.0089716796875,0.0000126953125,0];
semilogy(SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),R_MUMIMO(1,:)','b*--',...
         SNRMin+10*log10(coding_rate):SNRInc:SNRMax+10*log10(coding_rate),R_FSSMUMIMO(1,:)','rx--');
legend('MU-MIMO',...
       'FSS-MU-MIMO');
xlabel('Eb/No [dB]');
ylabel('BER');
axis([0,30,10^(-5),10^(0)]);
grid;