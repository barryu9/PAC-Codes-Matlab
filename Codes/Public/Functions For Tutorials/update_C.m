function C = update_C(pc_params, phi, C, u)
    N = pc_params.N;
    phi_mod_2 = mod(phi, 2);
    C(1, 1+phi_mod_2) = u;
    if (phi_mod_2 == 1) && (phi ~= N - 1)
        layer = pc_params.bit_layer_vec(phi + 1);
        for i_layer = 0:layer - 1
            index_1 = pc_params.lambda_offset(i_layer+1);
            index_2 = pc_params.lambda_offset(i_layer+2);
            for beta = index_1:2 * index_1 - 1
                C(beta + index_1, 2) = mod(C(beta, 1)+C(beta, 2), 2); %Left Column lazy copy
                C(beta + index_2, 2) = C(beta, 2);
    %             operation_count_C=operation_count_C+1;
    
            end
        end
        index_1 = pc_params.lambda_offset(layer+1);
        index_2 = pc_params.lambda_offset(layer+2);
        for beta = index_1:2 * index_1 - 1
            C(beta + index_1, 1) = mod(C(beta, 1)+C(beta, 2), 2); %Left Column lazy copy
            C(beta + index_2, 1) = C(beta, 2);
    %         operation_count_C=operation_count_C+1;
    
        end
    end
end