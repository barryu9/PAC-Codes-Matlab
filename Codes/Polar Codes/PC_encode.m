function x = pc_encode(pc_params,rp, d)
%encode Polar码编码
%将输入的源信息d编码成PAC码

    if (size(d,1) ~= pc_params.k)
        error('The length of the input d is not equal to k.')
    end
    batchsize=size(d,2);
    crc_parity_check=pc_params.crc_parity_check;
    info_bits_indices=rp.info_bits_indices;
    G_N=pc_params.G_N;
    N=pc_params.N;
    x=zeros(N,batchsize);

    for batch_index =1:batchsize
        % Rate Profile
        if (pc_params.crc_length > 0)
            info_with_crc = [d(:,batch_index); mod(crc_parity_check*d(:,batch_index), 2)];
        else
            info_with_crc = d(:,batch_index);
        end

        u = zeros(1, N);
        u(info_bits_indices) = info_with_crc;
        % Polar Encoding
        x(:,batch_index) = mod(u*G_N, 2)';
    end
end