function [g_crc] = get_crc_objective(crc_length)
switch crc_length
    case 4
        g_crc = [1, 0, 0, 1, 1];
    case 6
        g_crc = [1, 0, 0, 0, 0, 1, 1];
    case 8
        g_crc = [1, 1, 1, 1, 1, 1, 0, 0, 1];
    case 10
        g_crc = [1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1];
    case 12
        g_crc = [1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1];
    case 16
        g_crc = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1];
    case 24
        g_crc = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1];
    otherwise
        disp('Unsupported CRC length. Program terminates')
end