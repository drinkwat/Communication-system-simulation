% CDMA 系统仿真参数设置
num_users = 4; % 用户数量
chip_length = 100; % 码片长度
signal_length = 1000; % 信号长度
SNR_dB = 10; % 信噪比（以 dB 为单位）

% 生成伪随机码（Walsh 码，这里假设使用简单的正交码）
Walsh_code = cell(num_users, 1);
for i = 1:num_users
    Walsh_code{i} = randi([-1, 1], 1, chip_length);
end

% 生成用户发送的数据信号
data_signal = cell(num_users, 1);
for i = 1:num_users
    data = randi([0, 1], 1, signal_length);
    data_signal{i} = data;
end

% 对每个用户的数据信号进行扩频操作
spread_signal = cell(num_users, 1);
for i = 1:num_users
    spread_signal{i} = Walsh_code{i}.* data_signal{i};
end

% 叠加所有用户的扩频信号
combined_signal = zeros(1, chip_length * signal_length);
for i = 1:num_users
    combined_signal = combined_signal + spread_signal{i};
end

% 加入噪声（根据给定的 SNR）
SNR = 10^(SNR_dB/10);
noise_power = var(combined_signal) / SNR;
noise = sqrt(noise_power) * randn(1, chip_length * signal_length);
received_signal = combined_signal + noise;

% 接收端解扩
received_data = cell(num_users, 1);
for i = 1:num_users
    received_data{i} = received_signal.* Walsh_code{i};
end

% 判决检测（简单的阈值判决）
threshold = 0;
detected_data = cell(num_users, 1);
for i = 1:num_users
    detected_data{i} = received_data{i} > threshold;
end

% 计算误码率
total_errors = 0;
for i = 1:num_users
    errors = sum(abs(data_signal{i} - detected_data{i}));
    total_errors = total_errors + errors;
end
num_bits = sum(cellfun(@numel, data_signal));
BER = total_errors / num_bits;

disp(['误码率 (BER): ', num2str(BER)]);