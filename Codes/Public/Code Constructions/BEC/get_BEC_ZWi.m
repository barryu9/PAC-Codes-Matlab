function ZWi = get_BEC_ZWi(N, ZW)
ZWi = zeros(N, 1);
ZWi(1) = ZW;
for level = 1:log2(N)
    B = 2^level;
    for j = 1:B/2
        Z_tmp = ZWi(j);
        ZWi(j) = 2 * Z_tmp - Z_tmp^2; %use upper bound
        ZWi(B/2+j) = Z_tmp^2;
    end
end
ZWi = bitrevorder(ZWi); %perform bit-reversal order.