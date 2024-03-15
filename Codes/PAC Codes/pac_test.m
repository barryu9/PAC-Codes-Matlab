% g = [1,0,1,1,0,1,1]
% generate_CS([6,7,8,11,12,13,14,15,16],4)

clear
addpath(genpath('Codes/'))

N = 128;
k = 64;
crc_length = 0;
F_N=[1 0;1 1];
gen = [1 0 1 1 0 1,1];
pac_params = pac_init(N,k,crc_length,F_N,gen);
rp = rp_RM_Polar(N,k+crc_length,3);
snr_dB=3;
tic
for i = 1:1000
    u = double(rand(k,1)>0.5);
    x = pac_encode(pac_params,rp,u);
    sigma = 1/sqrt(2 * pac_params.R) * 10^(-snr_dB/20);
    bpsk = 1 - 2 * x;
    noise = randn(N, 1);
    y = bpsk + sigma * noise;
    llr = 2/sigma^2*y;
    d = pac_fano_decoder(pac_params,rp, llr, 1);
    if (sum(sum(d~=u))>0)
        disp("worse")
    end
end
toc