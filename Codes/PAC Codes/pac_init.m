function [pac_params] = pac_init(N,k,crc_length,F_N,gen)
%INIT_PAC_CODES 设置原始设定变量来初始化Polar码
%   输入：
%    - N 码长
%    - k 信息序列长度
%    - crc_length CRC编码长度
%    - F_N 极化核
%    - gen 卷积编码冲激序列

%  输出：
%   - pac_params 结构体，包含以下字段：
%       -- N 码长
%       -- n log2(N)
%       -- k 信息序列长度
%       -- R 码率 k/N
%       -- crc_length CRC编码长度
%       -- F_N 极化核
%       -- G_N 极化码生成矩阵
%       -- conv_depth 卷积深度
%       -- T   卷积预编码矩阵
%       -- lambda_offset
%       -- llr_layer_vec
%       -- bit_layer_vec
%       -- G_crc CRC生成矩阵
%       -- H_crc CRC校验矩阵

pac_params = struct;
pac_params.N = N;
pac_params.n = log2(N);
pac_params.k = k;
pac_params.R = k/N;
pac_params.crc_length = crc_length;
pac_params.F_N = F_N;
pac_params.G_N = get_G_N(F_N,N);
pac_params.lambda_offset = 2.^(0:pac_params.n);
pac_params.llr_layer_vec = get_llr_layer(N);
pac_params.bit_layer_vec = get_bit_layer(N);
pac_params.gen = gen;
pac_params.conv_depth = length(gen);

gen_zp = [gen, zeros(1, N-pac_params.conv_depth)];
pac_params.T = triu(toeplitz(gen_zp)); %upper-triangular Toeplitz matrix

if (crc_length > 0)
    g_crc = get_crc_objective(crc_length);
    % g_crc = [1, 0, 0, 1, 1]; 也可以自己设置g_crc;
    [G_crc, pac_params.H_crc] = crc_generator_matrix(g_crc, k);
    pac_params.crc_parity_check = G_crc(:, k+1:end)';
else
    pac_params.crc_parity_check = [];
    pac_params.H_crc = [];
end

end

