function P = update_P(pc_params, llr, phi, P, C)
N = pc_params.N;
n = log2(N);
layer = pc_params.llr_layer_vec(phi + 1);

switch phi %Decoding bits u_0 and u_N/2 needs channel LLR, so the decoding of them is separated from other bits.
    case 0
        index_1 = pc_params.lambda_offset(n);
        for beta = 0:index_1 - 1
            P(beta + index_1) = sign(llr(beta + 1)) * sign(llr(beta + index_1 + 1)) * min(abs(llr(beta + 1)), abs(llr(beta + index_1 + 1)));
%             operation_count_P=operation_count_P+1;
        end
        for i_layer = n - 2:-1:0
            index_1 = pc_params.lambda_offset(i_layer+1);
            index_2 = pc_params.lambda_offset(i_layer+2);
            for beta = 0:index_1 - 1
                P(beta + index_1) = sign(P(beta + index_2)) * ...
                    sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), abs(P(beta + index_1 + index_2)));
%                 operation_count_P=operation_count_P+1;
            end
        end
    case N / 2
        index_1 = pc_params.lambda_offset(n);
        for beta = 0:index_1 - 1
            x_tmp = C(beta + index_1, 1);
            P(beta + index_1) = (1 - 2 * x_tmp) * llr(beta + 1) + llr(beta + 1 + index_1);
%             operation_count_P=operation_count_P+1;

        end
        for i_layer = n - 2:-1:0
            index_1 = pc_params.lambda_offset(i_layer+1);
            index_2 = pc_params.lambda_offset(i_layer+2);
            for beta = 0:index_1 - 1
                P(beta + index_1) = sign(P(beta + index_2)) * ...
                    sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), abs(P(beta + index_1 + index_2)));
%                 operation_count_P=operation_count_P+1;

            end
        end
    otherwise
        index_1 = pc_params.lambda_offset(layer+1);
        index_2 = pc_params.lambda_offset(layer+2);
        for beta = 0:index_1 - 1
            P(beta + index_1) = (1 - 2 * C(beta + index_1, 1)) * P(beta + index_2) + ...
                P(beta + index_1 + index_2);
%             operation_count_P=operation_count_P+1;

        end
        for i_layer = layer - 1:-1:0
            index_1 = pc_params.lambda_offset(i_layer+1);
            index_2 = pc_params.lambda_offset(i_layer+2);
            for beta = 0:index_1 - 1
                P(beta + index_1) = sign(P(beta + index_2)) * ...
                    sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), ...
                    abs(P(beta + index_1 + index_2)));
%                 operation_count_P=operation_count_P+1;

            end
        end
end
end
