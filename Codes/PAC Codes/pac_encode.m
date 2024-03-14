function x = pac_encode(pac_params,rp, d)
    %encode PAC编码
    %将输入的源信息d编码成PAC码
    if (size(d,1) ~= pac_params.k)
        error('The length of the input d is not equal to k.')
    end
    batchsize=size(d,2);
    crc_parity_check=pac_params.crc_parity_check;
    info_bits_indices=rp.info_bits_indices;
    G_N=pac_params.G_N;
    N=pac_params.N;
    x=zeros(N,batchsize);
    T=pac_params.T;
    for batch_index =1:batchsize
    
        % Rate Profile
        if (pac_params.crc_length > 0)
            info_with_crc = [d(:,batch_index); mod(crc_parity_check*d(:,batch_index), 2)];
        else
            info_with_crc = d(:,batch_index);
        end
        v = zeros(1, N);
        v(info_bits_indices) = info_with_crc;
        % convolutional encoder
        u = mod(v*T, 2);
        % Polar Encoding
        x(:,batch_index) = mod(u*G_N, 2)';
    end
end