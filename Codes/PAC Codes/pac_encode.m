function x = pac_encode(pac_params,rp, d)
    %encode PAC编码
    %将输入的源信息d编码成PAC码
    if (length(d) ~= pac_params.k)
        error('The length of the input d is not equal to k.')
    end
    % Rate Profile
    if (pac_params.crc_length > 0)
        info_with_crc = [d; mod(pac_params.crc_parity_check*d, 2)];
    else
        info_with_crc = d;
    end
    v = zeros(1, pac_params.N);
    v(rp.info_bits_indices) = info_with_crc;
    % convolutional encoder
    u = mod(v*pac_params.T, 2);
    % Polar Encoding
    x = mod(u*pac_params.G_N, 2)';
end