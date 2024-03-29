function rp = rp_modify(rp,new_info_bits_indices,indices_to_be_replaced)
    
    % 本函数不会对修改后的内容作校验，有可能修改前后 k 的大小会发生变化。

    info_bits_mask = rp.info_bits_mask;
    info_bits_mask(new_info_bits_indices) = 1;
    info_bits_mask(indices_to_be_replaced) = 0;
    
    info_bits_indices = find(info_bits_mask==1);
    frozen_bits_mask = 1 - info_bits_mask;
    frozen_bits_indices = find(frozen_bits_mask == 1);
    
    rp.info_bits_indices=info_bits_indices;
    rp.info_bits_mask=info_bits_mask;
    rp.frozen_bits_mask=frozen_bits_mask;
    rp.frozen_bits_indices=frozen_bits_indices;

end

