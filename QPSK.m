clc;
clear;
close all;

f = 2000;
fs = 1000 * 16;
snr = 20;%信噪比

%生成序列
I = 2 * randi([0 1],1,400) - 1;
Q = 2 * randi([0 1],1,400) - 1;

%上采样
I_up = reshape([I', zeros(length(I), 15)]', 1, 16 * length(I));
Q_up = reshape([Q', zeros(length(Q), 15)]', 1, 16 * length(Q));
%成型滤波
rcos_fir = rcosdesign(0.5, 6, 16);
rcos_I_up = conv(I_up, rcos_fir);
rcos_Q_up = conv(Q_up, rcos_fir);
filter_delay1 = (length(rcos_fir) - 1) / 2;
%调制
send = rcos_I_up.*cos(2*pi*f.*(1:length(rcos_I_up)) / fs) - rcos_Q_up.*sin(2*pi*f.*(1:length(rcos_Q_up)) / fs);
figure(1);
subplot(2,1,1);
plot(send(filter_delay1 + 1: filter_delay1 + 1 + 160));
title("发送的时域信号");
subplot(2,1,2);
plot(abs(fft(send, length(send))));
title("发送的频域信号");


%信道
send_noise = awgn(send, snr);
%

%相干解调
rec_I = send_noise.*cos(2*pi*f.*(1:length(send)) / fs);
rec_Q = send_noise.*sin(2*pi*f.*(1:length(send)) / fs);
%低通滤波
fir_lp =fir1(128,0.2); %截止频率为0.2*(fs/2) 使用汉明窗设计一个128阶低通带线性相位的FIR滤波器。
rec_I_lp = conv(fir_lp,rec_I);
rec_Q_lp = conv(fir_lp,rec_Q);
filter_delay2 = (length(fir_lp)-1)/2;  %低通滤波器的延迟时长
%匹配滤波
rcos_rec_I_lp = conv(rec_I_lp, rcos_fir);
rcos_rec_Q_lp = conv(rec_Q_lp, rcos_fir);
filter_delay3 = (length(rcos_fir) - 1) / 2;
%滤波器的延迟
total_delay = filter_delay1 + filter_delay2 + filter_delay3;
com_I = (rcos_rec_I_lp(total_delay + 1: 16 :end - total_delay));
com_Q = -(rcos_rec_Q_lp(total_delay + 1: 16 :end - total_delay));

IQ = complex(com_I + 1j .* com_Q);

scatterplot(IQ);title("QPSK星座图");
%判决
judge_I = sign(com_I);
judge_Q = sign(com_Q);

wrong = sum(~(judge_I==I)) + sum(~(judge_Q==Q));
fprintf("误码：%d\n", wrong);

% ebn0 = -6:8;
% ber = berawgn(ebn0,'psk',4,'nondiff');

%send = awgn(send, 8);
% rec = send .* sin(2*pi*f.*t / fs);
% rec_2 = sign(rec);
% 
% figure(2);
% plot(rec_2);
