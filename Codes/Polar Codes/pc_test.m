% g = [1,0,1,1,0,1,1]
% generate_CS([6,7,8,11,12,13,14,15,16],4)

clear
addpath(genpath('Codes/'))

N = 512;
k = 256;
crc_length = 4;
F_N=[1 0;1 1];

pc_params = pc_init(N,k,crc_length,F_N);
rp = GA_rate_profiling(N,k+crc_length,3);
snr_dB=3;
tic
for i = 1:1000
    u = double(rand(k,1)>0.5);
    x = pc_encode(pc_params,rp,u);
    sigma = 1/sqrt(2 * pc_params.R) * 10^(-snr_dB/20);
    bpsk = 1 - 2 * x;
    noise = randn(N, 1);
    y = bpsk + sigma * noise;
    llr = 2/sigma^2*y;
    d = pc_SCL_decoder(pc_params,rp, llr, 128);
    if (sum(sum(d~=u))>0)
        display("worse")
    end
end
toc