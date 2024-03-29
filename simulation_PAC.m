clc
clear
addpath(genpath('..'))

N = 128;
k = 64;
crc_length = 0;

F_N=[1 0;1 1];
gen = [1 0 1 1 0 1 1];

pac_params = pac_init(N,k,crc_length,F_N,gen);

dsnr_dB = 2.5;
rp = rp_RM_Polar(N,k+crc_length,dsnr_dB);
% rp = rp_modify(rp,[91,93,102,103,106,143,150],[173,179,203,213,226,233,241]);

snr_dB = [1.5];
L=32; 

min_block_num = 100000;
max_block_num = 100000;
max_error_num = 100000;


frame_errors_count=zeros(1,length(snr_dB));
bit_errors_count=zeros(1,length(snr_dB));
n_iter = zeros(1,length(snr_dB));

update_frequency = 1; %命令行输出的频率
batchsize = 100; % 一次仿真多少个batch




for i=1:length(snr_dB)
    sigma = 1/sqrt(2 * pac_params.R) * 10^(-snr_dB(i)/20);
    elapsetime_total=0;
    fprintf("Now Running SNR(dB)=%f\n",snr_dB(i))
    for ii = 1:max_block_num/batchsize
        tic;

        % 自适应仿真次数，优先保证达到最小仿真次数(min_iterations)
        % 如果错误个数大于max_iterations则直接结束

        if ii < min_block_num/batchsize
        elseif frame_errors_count(i) < max_error_num
        else
            continue;
        end


        d = double(rand(k,batchsize)>0.5);
        x = pac_encode(pac_params,rp,d);
        bpsk = 1 - 2 * x;
        noise = randn(N, batchsize);
        y = bpsk + sigma * noise;
        llr = 2/sigma^2*y;
        d_esti = pac_SCL_decoder(pac_params,rp,llr,L);
        errs=sum(d~=d_esti);
        n_iter(i)=n_iter(i)+batchsize;
        frame_errors_count(i)=frame_errors_count(i)+sum(errs>0);
        bit_errors_count(i)=bit_errors_count(i)+sum(errs);

        update_mod = mod(ii,update_frequency);
        elapsetime_total=elapsetime_total+toc;
        elapsetime_average = elapsetime_total/n_iter(i);

        if(update_mod==0)
           fprintf("%.2fdB@%i, Block Error(s):%i, BLER=%.2e; Bit Error(s):%i, BER=%.2e; %.2f bk/s, %s remaining\n",...
              snr_dB(i),ii*batchsize,frame_errors_count(i),frame_errors_count(i)/n_iter(i),bit_errors_count(i),bit_errors_count(i)/(n_iter(i)*k),...
              1/elapsetime_average, string(seconds(elapsetime_average*batchsize*(max_block_num/batchsize-ii)),"hh:mm:ss"))
        end
    end

end


FER=frame_errors_count./n_iter;
BER=bit_errors_count./(n_iter.*k);

save(['results\PAC_',datestr(datetime('now'),'yyyy-mm-dd-HH-MM'),'.mat'])

figure;
semilogy(snr_dB,FER,'-o','LineWidth',1);
grid on;



