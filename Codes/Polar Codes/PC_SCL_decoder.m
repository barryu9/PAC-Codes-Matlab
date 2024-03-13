function d_esti = PC_SCL_decoder(pc_params,rp, llr, L)
%LLR-based SCL deocoder, a single function, no other sub-functions.
%Frequently calling sub-functions will derease the efficiency of MATLAB
%codes.
%const
local N;
local k;
local lambda_offset;
local bit_layer_vec;
local llr_layer_vec;
local phi;
local C;
local P;

N = pc_params.N;
k = pc_params.k;
frozen_bits_mask = rp.frozen_bits_mask;
lambda_offset = pc_params.lambda_offset;
llr_layer_vec = pc_params.llr_layer_vec;
bit_layer_vec = pc_params.bit_layer_vec;
H_crc = pc_params.H_crc;
crc_length = pc_params.crc_length;
%memory declared
%If you can understand lazy copy and you just start learning polar codes
%for just fews days, you are very clever,
P = zeros(N - 1, L); %Channel llr is public-used, so N - 1 is enough.
C = zeros(N - 1, 2 * L); %I do not esitimate (x1, x2, ... , xN), so N - 1 is enough.
d = zeros(k, L); %unfrozen bits that polar codes carry, including crc bits.
PM = zeros(L, 1); %Path metrics
activepath = zeros(L, 1); %Indicate if the path is active. '1'→active; '0' otherwise.
cnt_u = 1; %information bit counter
 
%initialize
activepath(1) = 1;
%decoding starts
%default: in the case of path clone, the origianl path always corresponds to bit 0, while the new path bit 1.
for phi = 0:N - 1
    for l_index = 1:L
        if activepath(l_index) == 0
            continue;
        end
        P(:, l_index) = update_P_test(N,lambda_offset,llr_layer_vec,llr, phi, P(:, l_index), C(:, 2*l_index-1:2*l_index));
    end
    if frozen_bits_mask(phi + 1) == 0 %if now we decode an unfrozen bit
        PM_pair = realmax * ones(2, L);
        for l_index = 1:L
            if activepath(l_index) == 0
                continue;
            end
            PM_pair(1, l_index) = calc_PM(PM(l_index), P(1, l_index), 0);
            PM_pair(2, l_index) = calc_PM(PM(l_index), P(1, l_index), 1);
        end
        middle = min(2*sum(activepath), L);
        PM_sort = sort(PM_pair(:));
        PM_cv = PM_sort(middle);
        compare = PM_pair <= PM_cv;
        kill_index = zeros(L, 1); %to record the index of the path that is killed
        kill_cnt = 0; %the total number of killed path
        %the above two variables consist of a stack
        for i = 1:L
            if (compare(1, i) == 0) && (compare(2, i) == 0) %which indicates that this path should be killed
                activepath(i) = 0;
                kill_cnt = kill_cnt + 1; %push stack
                kill_index(kill_cnt) = i;
            end
        end
        for l_index = 1:L
            if activepath(l_index) == 0
                continue;
            end
            path_state = compare(1, l_index) * 2 + compare(2, l_index);
            switch path_state %path_state can equal to 0, but in this case we do no operation.
                case 1 % PM of the second row is lower
                    d(cnt_u, l_index) = 1;
                    C(:, 2*l_index-1:2*l_index) = update_C_test(N,lambda_offset,bit_layer_vec, phi, C(:, 2*l_index-1:2*l_index), 1);
                    PM(l_index) = PM_pair(2, l_index);
                case 2 % PM of the first row is lower
                    d(cnt_u, l_index) = 0;
                    C(:, 2*l_index-1:2*l_index) = update_C_test(N,lambda_offset,bit_layer_vec, phi, C(:, 2*l_index-1:2*l_index), 0);
                    PM(l_index) = PM_pair(1, l_index);
                case 3 %
                    index = kill_index(kill_cnt);
                    kill_cnt = kill_cnt - 1; %pop stack
                    activepath(index) = 1;
                    %lazy copy
                    C(:, 2*index-1:2*index) = C(:, 2*l_index-1:2*l_index);
                    P(:, index) = P(:, l_index);
                    d(:, index) = d(:, l_index);
                    d(cnt_u, l_index) = 0;
                    d(cnt_u, index) = 1;
                    C(:, 2*l_index-1:2*l_index) = update_C_test(N,lambda_offset,bit_layer_vec, phi, C(:, 2*l_index-1:2*l_index), 0);
                    C(:, 2*index-1:2*index) = update_C_test(N,lambda_offset,bit_layer_vec, phi, C(:, 2*index-1:2*index), 1);
                    PM(l_index) = PM_pair(1, l_index);
                    PM(index) = PM_pair(2, l_index);
            end
        end
        cnt_u = cnt_u + 1;
    else %frozen bit operation
        for l_index = 1:L
            if activepath(l_index) == 0
                continue;
            end
            PM(l_index) = calc_PM(PM(l_index), P(1, l_index), 0);
            C(:, 2*l_index-1:2*l_index) = update_C_test(N,lambda_offset,bit_layer_vec, phi, C(:, 2*l_index-1:2*l_index), 0);

        end
    end


end
%path selection.
activepath = logical(activepath);
PM_active = PM(activepath);
d_active = d(:, activepath);
[~, path_ordered] = sort(PM_active);
if (crc_length > 0)
    for l_index = 1:length(PM_active)
        path_num = path_ordered(l_index);
        info_with_crc = d(:, path_num);
        err = sum(mod(H_crc*info_with_crc, 2));
        if err == 0
            d_esti = d_active(1:end-crc_length, path_num);
            break;
        else
            if l_index == length(PM_active)
                d_esti = d_active(1:end-crc_length, path_ordered(1));
            end
        end
    end
else
    d_esti = d_active(1:end-crc_length, path_ordered(1));
end
end

function P = update_P_test(N,lambda_offset,llr_layer_vec, llr, phi, P, C)

n = log2(N);
layer = llr_layer_vec(phi + 1);

switch phi %Decoding bits u_0 and u_N/2 needs channel LLR, so the decoding of them is separated from other bits.
    case 0
        index_1 = lambda_offset(n);
        for beta = 0:index_1 - 1
            P(beta + index_1) = sign(llr(beta + 1)) * sign(llr(beta + index_1 + 1)) * min(abs(llr(beta + 1)), abs(llr(beta + index_1 + 1)));
%             operation_count_P=operation_count_P+1;
        end
        for i_layer = n - 2:-1:0
            index_1 = lambda_offset(i_layer+1);
            index_2 = lambda_offset(i_layer+2);
            for beta = 0:index_1 - 1
                P(beta + index_1) = sign(P(beta + index_2)) * ...
                    sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), abs(P(beta + index_1 + index_2)));
%                 operation_count_P=operation_count_P+1;
            end
        end
    case N / 2
        index_1 = lambda_offset(n);
        for beta = 0:index_1 - 1
            x_tmp = C(beta + index_1, 1);
            P(beta + index_1) = (1 - 2 * x_tmp) * llr(beta + 1) + llr(beta + 1 + index_1);
%             operation_count_P=operation_count_P+1;

        end
        for i_layer = n - 2:-1:0
            index_1 = lambda_offset(i_layer+1);
            index_2 = lambda_offset(i_layer+2);
            for beta = 0:index_1 - 1
                P(beta + index_1) = sign(P(beta + index_2)) * ...
                    sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), abs(P(beta + index_1 + index_2)));
%                 operation_count_P=operation_count_P+1;

            end
        end
    otherwise
        index_1 = lambda_offset(layer+1);
        index_2 = lambda_offset(layer+2);
        for beta = 0:index_1 - 1
            P(beta + index_1) = (1 - 2 * C(beta + index_1, 1)) * P(beta + index_2) + ...
                P(beta + index_1 + index_2);
%             operation_count_P=operation_count_P+1;

        end
        for i_layer = layer - 1:-1:0
            index_1 = lambda_offset(i_layer+1);
            index_2 = lambda_offset(i_layer+2);
            for beta = 0:index_1 - 1
                P(beta + index_1) = sign(P(beta + index_2)) * ...
                    sign(P(beta + index_1 + index_2)) * min(abs(P(beta + index_2)), ...
                    abs(P(beta + index_1 + index_2)));
%                 operation_count_P=operation_count_P+1;

            end
        end
end
end

function C = update_C_test(N,lambda_offset,bit_layer_vec, phi, C, u)
local N;
local k;
local lambda_offset;
local bit_layer_vec;
local llr_layer_vec;
local phi;
local C;
local P;

phi_mod_2 = mod(phi, 2);
C(1, 1+phi_mod_2) = u;
if (phi_mod_2 == 1) && (phi ~= N - 1)
    layer = bit_layer_vec(phi + 1);
    for i_layer = 0:layer - 1
        index_1 = lambda_offset(i_layer+1);
        index_2 = lambda_offset(i_layer+2);
        for beta = index_1:2 * index_1 - 1
            C(beta + index_1, 2) = mod(C(beta, 1)+C(beta, 2), 2); %Left Column lazy copy
            C(beta + index_2, 2) = C(beta, 2);
%             operation_count_C=operation_count_C+1;

        end
    end
    index_1 = lambda_offset(layer+1);
    index_2 = lambda_offset(layer+2);
    for beta = index_1:2 * index_1 - 1
        C(beta + index_1, 1) = mod(C(beta, 1)+C(beta, 2), 2); %Left Column lazy copy
        C(beta + index_2, 1) = C(beta, 2);
%         operation_count_C=operation_count_C+1;

    end
end
end