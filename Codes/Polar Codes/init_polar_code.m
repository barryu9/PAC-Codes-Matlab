function [pc_params] = init_polar_code(N,k,crc_length,F_N)
%INIT_POLAR_CODES 设置原始设定变量来初始化Polar码
%   输入：
%    - N 码长
%    - k 信息序列长度
%    - crc_length CRC编码长度
%    - F_N 极化核

%  输出：
%   - pc_paras 结构体，包含以下字段：
%       -- N 码长
%       -- n log2(N)
%       -- k 信息序列长度
%       -- R 码率 k/N
%       -- crc_length CRC编码长度
%       -- F_N 极化核
%       -- G_N 极化码生成矩阵
%       -- lambda_offset
%       -- llr_layer_vec
%       -- bit_layer_vec
%       -- G_crc CRC生成矩阵
%       -- H_crc CRC校验矩阵

pc_params = struct;
pc_params.N = N;
pc_params.k = k;
pc_params.R = k/N;
pc_params.crc_length = crc_length;
pc_params.F_N = F_N;

pc_params.n = log2(N);
pc_params.G_N = get_G_N(F_N,N);
pc_params.lambda_offset = 2.^(0:pc_params.n);
pc_params.llr_layer_vec = get_llr_layer(N);
pc_params.bit_layer_vec = get_bit_layer(N);

if (crc_length > 0)
    g_crc = get_crc_objective(crc_length);
    % g_crc = [1, 0, 0, 1, 1]; 也可以自己设置g_crc;
    [G_crc, pc_params.H_crc] = crc_generator_matrix(g_crc, k);
    pc_params.crc_parity_check = G_crc(:, k+1:end)';
else
    pc_params.crc_parity_check = [];
    pc_params.H_crc = [];
end

end

