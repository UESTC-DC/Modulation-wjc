clc;
close all;
clear;

%%%%%%%%%%%参数设定%%%%%%%%%%%%%%%%%%%%

bit_rate = 1000;            %比特率
% symbol_rate = 1000;         %符号率
sample = 16;                %每个符号的采样点数
fc = 2000;                  %载波频率
fs = bit_rate*sample;       %采样频率=比特率*每个符号的采样点数
source_number = 352;        %发送信号的长度
rollof_factor = 0.5;        %滚降因子,可调整

%%%%%%%%%%%%%%%%信源%%%%%%%%%%%%%%%%%%

%随机信号
source = randi([0 1],1,source_number);
%给出标志性帧头，方便调试
frame_pre = ones(1,32); %用于捕获和同步
frame_begin = [0 1 1 1 1 1 1 0]; %帧开始的标志
frame_head = [frame_pre frame_begin];  %帧数据
frame_end = [0 1 1 1 1 1 1 0];   %帧结束的标志
%组帧
frame_msg = [frame_head source frame_end];

%%%%%%%%%%%%发射机%%%%%%%%%%%%%%%%%%%%

%双极性转换
bipolar_msg_source = 2*frame_msg-1;   %双极性信号
% 波形观察
figure(1);
subplot(211);plot(bipolar_msg_source);
title('时域波形');
subplot(212);plot(abs(fft(bipolar_msg_source)));
title('频域波形');

%上采样
bipolar_msg_source_temp = [bipolar_msg_source',zeros(size(bipolar_msg_source,2),sample-1)];
length_x = size(bipolar_msg_source_temp,1);
length_y = size(bipolar_msg_source_temp,2);
up16_bipolar_msg_source = reshape(bipolar_msg_source_temp',1,length_x*length_y);

%波形观察
figure(2);
subplot(211);plot(up16_bipolar_msg_source);
title('上采样时域波形');
subplot(212);plot(abs(fft(up16_bipolar_msg_source)));
title('上采样频域波形');

%%%滤波器
%滚降滤波器
rcos_fir = rcosdesign(rollof_factor,6,sample);
fvtool(rcos_fir,'Analysis','impulse');    %将脉冲响应可视化
%采用滚降滤波器进行滤波
% rcos_msg_source = filter(rcos_fir,1,up16_bipolar_msg_source);
rcos_msg_source = conv(up16_bipolar_msg_source,rcos_fir);
filter_delay1 = (length(rcos_fir)-1)/2;  %滚降滤波器的延迟时长
% rcos_msg_source_TEMP = zeros(1,length(rcos_msg_source));
% Z = [rcos_msg_source',rcos_msg_source_TEMP'];
% Z=single(Z);
% charSymbolrate=num2str(bit_rate/1000);
% charFreq_DEV=num2str(fs/1000);
% csvName = ['BPSK_SymbolRate',charSymbolrate,'kbps_FreqDev',charFreq_DEV, 'kHz.csv'];
% csvwrite(csvName,Z,0,0);
%波形观察
figure(3);
subplot(211);plot(rcos_msg_source(1:1024));
title('通过成型滤波器的时域波形');
subplot(212);plot(abs(fft(rcos_msg_source)));
title('通过成型滤波器的频域波形');

%%%%%%%%%%%%调制器%%%%%%%%%%%%%%%%%%%%
%%%载波发送
time = 1:length(rcos_msg_source);
rcos_msg_source_carrier = rcos_msg_source.*cos(2*pi*fc.*time/fs);
%波形观察
figure(4);
subplot(211);plot(rcos_msg_source_carrier(1:1024));
title('载波调制时域波形');
subplot(212);plot(abs(fft(rcos_msg_source_carrier)));
title('载波调制频域波形');

%%%%%%%%%%%%信道%%%%%%%%%%%%%%%%%%%%%%
%设置信噪比
ebn0 = -6:8;
snr = ebn0 - 10*log10(0.5*16);
for i = 1:length(snr)
    %线性高斯白噪声信道
    rcos_msg_source_carrier_addnoise = awgn(rcos_msg_source_carrier,snr(i),'measured');
%     %波形观察
%     figure(11);
%     plot(rcos_msg_source_carrier_addnoise);
%     title('时域波形');
%     figure(12);
%     plot(abs(fft(rcos_msg_source_carrier_addnoise)));
%     title('频域波形');
    
    %%%%%%%%%%%%接收机%%%%%%%%%%%%%%%%%%%%
    %%%%%%载波恢复
    %%%相干解调
    rcos_msg_source_addnoise =rcos_msg_source_carrier_addnoise.*cos(2*pi*fc.*time/fs);

%     %波形观察
%     figure(13);
%     plot(rcos_msg_source_addnoise);
%     title('时域波形');
%     figure(14);
%     plot(abs(fft(rcos_msg_source_addnoise)));
%     title('频域波形');


    %%%%%%%滤波
    %%%%低通滤波
    fir_lp =fir1(128,0.2); %截止频率为0.2*(fs/2) 使用汉明窗设计一个128阶低通带线性相位的FIR滤波器。
    rcos_msg_source_1p = conv(fir_lp,rcos_msg_source_addnoise);
    %延迟64个采样点输出
    filter_delay2 = (length(fir_lp)-1)/2;  %低通滤波器的延迟时长

%     %波形观察
%     figure(15);
%     plot(rcos_msg_source_lp);
%     title('时域波形');
%     figure(16);
%     plot(abs(fft(rcos_msg_source_lp)));
%     title('频域波形');

    %%%%%%匹配滤波
    %生成匹配滤波器
    rcos_fir = rcosdesign(rollof_factor,6,sample);
    %滤波
    rcos_msg_source_MF = conv(rcos_fir,rcos_msg_source_1p);
    filter_delay3 = (length(rcos_fir)-1)/2;  %滚降滤波器的延迟时长
%     %波形观察
%     figure(17);
%     plot(rcos_msg_source_MF);
%     title('时域波形');
%     figure(18);
%     plot(abs(fft(rcos_msg_source_MF)));
%     title('频域波形');

    %%%%%最佳采样
    %%%选取最佳采样点
    decision_site = filter_delay1+filter_delay2+filter_delay3; %(96+128+96)/2 =160 三个滤波器的延迟 96 128 96


    %每个符号选取一个点作为判决
    rcos_msg_source_MF_option = rcos_msg_source_MF(decision_site+1:sample:end-decision_site);
    %涉及到三个滤波器，固含有滤波器延迟累加


    %%判决
    msg_source_MF_option_sign= sign(rcos_msg_source_MF_option);

%     %波形观察
%     figure(13);
%     plot(msg_source_MF_option_sign,'-*');
%     title('判决结果');
%     
%     eyediagram(rcos_msg_source,sample);
%     title('发射端眼图');
%     eyediagram(rcos_msg_source_MF,sample);
%     title('接收端眼图');
%     
%     scatterplot(rcos_msg_source(48+1:16:end-48));
%     title('BPSK星座图');


%%%%%%%%%%%%%%%%%   信宿    %%%%%%%%%%%%%%%%%%%%
%%%误码率性能比对
%[err_number,bit_err_ratio]=biterr(x,y)
[err_number(i),bit_err_ratio(i)]=biterr(frame_msg(1:length(rcos_msg_source_MF_option)),(msg_source_MF_option_sign+1)/2);

end 
%%%%%%%%%%%%%%%%%   仿真结果    %%%%%%%%%%%%%%%%%%%%
ber = berawgn(ebn0,'psk',2,'nondiff');
figure(5);
semilogy(ebn0,bit_err_ratio,'-*',ebn0,ber,'-+');
xlabel('比特信噪比');
ylabel('误码率');
title('不同信噪比下误码率仿真曲线');
legend('实验曲线','理论曲线');
grid on;