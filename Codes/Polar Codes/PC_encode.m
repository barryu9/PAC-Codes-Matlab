function x = pc_encode(pc_params,rp, d)
%encode Polar码编码
%将输入的源信息d编码成PAC码
    if (length(d) ~= pc_params.k)
        error('The length of the input d is not equal to k.')
    end
    % Rate Profile
    if (pc_params.crc_length > 0)
        info_with_crc = [d; mod(pc_params.crc_parity_check*d, 2)];
    else
        info_with_crc = d;
    end
    u = zeros(1, pc_params.N);
    u(rp.info_bits_indices) = info_with_crc;
    % Polar Encoding
    x = mod(u*pc_params.G_N, 2)';
end