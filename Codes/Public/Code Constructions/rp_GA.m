function rp = rp_GA(N,k,dsnr_dB)
    %UNTITLED 此处显示有关此函数的摘要
    %   此处显示详细说明
    
    rp = struct;
    R = k/N;
    
    sigma = 1 / sqrt(2*R) * 10^(-dsnr_dB / 20);
    [channel_params, ~] = GA(sigma, N);
    [~, channels_ordered] = sort(channel_params, 'descend');
    
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

