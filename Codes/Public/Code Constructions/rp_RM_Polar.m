function rp = rp_RM_Polar(N,k,dsnr_dB)
    %UNTITLED 此处显示有关此函数的摘要
    %   此处显示详细说明
    
    rp = struct;
    R = k/N;
    channel_indices = (0:N - 1)';
    bitStr = dec2bin(channel_indices);
    bit = abs(bitStr) - 48;
    RM_score = sum(bit, 2);
    
    sigma = 1 / sqrt(2*R) * 10^(-dsnr_dB / 20);
    [channel_params, ~] = GA(sigma, N);
    [RM_score_ordered, ] = sort(RM_score, 'descend');    % 先将RM得分排序
    RM_score_k = RM_score_ordered(k); %找出第k位的RM得分
    channel_params_temp = channel_params;
    channel_params_temp(RM_score<RM_score_k) = -1; 
    %将所有RM得分小于RM_score_k的信道得分设为-1，也就是说排序永远在最下面，直接淘汰！
    channel_params_temp(RM_score>RM_score_k) = max(channel_params)+1; 
    %将所有RM得分大于RM_score_k的信道得分设为最大值+1，也就是说排序永远在最上面，直接晋级！

    [~, channels_ordered] = sort(channel_params_temp, 'descend');
    %这里实际上只对RM得分等于RM_score_k的信道进行排序
    
    info_bits_indices = sort(channels_ordered(1:k), 'ascend');
    info_bits_mask = zeros(1,N);
    info_bits_mask(info_bits_indices) = 1;
    frozen_bits_mask = 1 - info_bits_mask;
    frozen_bits_indices = find(frozen_bits_mask == 1);
    
    rp.info_bits_indices=info_bits_indices;
    rp.info_bits_mask=info_bits_mask;
    rp.frozen_bits_mask=frozen_bits_mask;
    rp.frozen_bits_indices=frozen_bits_indices;
    rp.channel_params = channel_params;
    rp.pe = 0.5-0.5*erf(sqrt(channel_params)/2);
    
end

