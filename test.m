% g = [1,0,1,1,0,1,1]
% generate_CS([6,7,8,11,12,13,14,15,16],4)

clear
addpath(genpath('Codes/'))

N = 512;
k = 256;
crc_length = 0;
F_N=[1 0;1 1];

pc_params = init_polar_code(N,k,crc_length,F_N);
rp = GA_rate_profiling(N,k,3);
snr_dB=100;
tic
for i = 1:1000
    u = double(rand(k,1)>0.5);
    x = PC_encode(pc_params,rp,u);
    sigma = 1/sqrt(2 * pc_params.R) * 10^(-snr_dB/20);
    bpsk = 1 - 2 * x;
    noise = randn(N, 1);
    y = bpsk + sigma * noise;
    llr = 2/sigma^2*y;
    d = PC_SCL_decoder(pc_params,rp, llr, 32);
%     d = PC_SCL_decoder_test(N,k,rp.frozen_bits_mask,crc_length,pc_params.H_crc,pc_params.lambda_offset,pc_params.llr_layer_vec,pc_params.bit_layer_vec,llr, 32);
    if d~=u
        printf("worse")
    end
end
toc