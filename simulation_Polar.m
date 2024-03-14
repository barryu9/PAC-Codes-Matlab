clc
clear
addpath(genpath('..'))

N = 256;
k = 128;
crc_length = 0;

F_N=[1 0;1 1];
pac_params = pc_init(N,k,crc_length,F_N);

dsnr_dB = 2.5;

rp = rp_GA(N,k+crc_length,dsnr_dB);
% rp = rp_modify(rp,[91,93,102,103,106,143,150],[173,179,203,213,226,233,241]);

snr_dB = [2.5];
L=32;

min_iterations = 10000;
max_iterations = 10000;
max_error_num = 10000;


frame_errors_count=zeros(1,length(snr_dB));
bit_errors_count=zeros(1,length(snr_dB));
n_iter = zeros(1,length(snr_dB));

update_frequency = 100; %命令行输出的频率
elapsetime_filter = zeros(1,update_frequency);


for i=1:length(snr_dB)
    sigma = 1/sqrt(2 * pac_params.R) * 10^(-snr_dB(i)/20);
    fprintf("Now Running SNR(dB)=%f\n",snr_dB(i))
    for ii = 1:max_iterations
        tic;

        % 自适应仿真次数，优先保证达到最小仿真次数(min_iterations)
        % 如果错误个数大于max_iterations则直接结束

        if ii < min_iterations
        elseif frame_errors_count(i) < max_error_num
        else
            continue;
        end


        u = double(rand(k,1)>0.5);
        x = pc_encode(pac_params,rp,u);
        bpsk = 1 - 2 * x;
        noise = randn(N, 1);
        y = bpsk + sigma * noise;
        llr = 2/sigma^2*y;
        d = pc_SCL_decoder(pac_params,rp,llr,L);
        errs=sum(sum(u~=d));
        n_iter(i)=ii;
        if(errs>0)
            frame_errors_count(i)=frame_errors_count(i)+1;
            bit_errors_count(i)=bit_errors_count(i)+errs;
        end

        update_mod = mod(ii,update_frequency);
        elapsetime_filter(update_mod+1)=toc;
        elapsetime_average = mean(elapsetime_filter);

        if(update_mod==0)
           fprintf("%.2fdB@%i, Block Error(s):%i, BLER=%.2e; Bit Error(s):%i, BER=%.2e; %.2f it/s, %s remaining\n",...
              snr_dB(i),ii,frame_errors_count(i),frame_errors_count(i)/ii,bit_errors_count(i),bit_errors_count(i)/ii,...
              1/elapsetime_average, string(seconds(elapsetime_average*(max_iterations-ii)),"hh:mm:ss"))
        end
    end

end


FER=frame_errors_count./n_iter;
BER=bit_errors_count./(n_iter.*k);

save(['results\PAC_',datestr(datetime('now'),'yyyy-mm-dd-HH-MM'),'.mat'])

figure;
semilogy(snr_dB,FER,'-o','LineWidth',1);
grid on;



