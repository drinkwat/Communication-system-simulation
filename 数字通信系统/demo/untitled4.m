function varargout = untitled4(varargin)
% UNTITLED4 MATLAB code for untitled4.fig
%      UNTITLED4, by itself, creates a new UNTITLED4 or raises the existing
%      singleton*.
%
%      H = UNTITLED4 returns the handle to a new UNTITLED4 or the handle to
%      the existing singleton*.
%
%      UNTITLED4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED4.M with the given input arguments.
%
%      UNTITLED4('Property','Value',...) creates a new UNTITLED4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled4

% Last Modified by GUIDE v2.5 02-Jan-2025 12:33:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled4_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled4_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before untitled4 is made visible.
function untitled4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled4 (see VARARGIN)

% Choose default command line output for untitled4
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

 % 基本参数定义
    handles.n = 15; % 产生码元数    
    handles.k = 7;
    handles.L = 100; % 每码元复制L次,每个码元采样次数
    handles.Ts = 0.001; % 每个码元的宽度
    handles.Rb = 1 / handles.Ts; % 码元速率1Kbps
    handles.dt = handles.Ts / handles.L; % 采样间隔
    handles.TotalT = handles.n * handles.Ts; % 总时间
    handles.t = 0:handles.dt:handles.TotalT - handles.dt; % 时间
    handles.Fs = 1 / handles.dt; % 采样频率

    % 校验矩阵
    handles.H = [
        1 1 0 1 0 0 0 1 0 0 0 0 0 0 0;
        0 1 1 0 1 0 0 0 1 0 0 0 0 0 0; 
        0 0 1 1 0 1 0 0 0 1 0 0 0 0 0; 
        0 0 0 1 1 0 1 0 0 0 1 0 0 0 0; 
        1 1 0 1 1 1 0 0 0 0 0 1 0 0 0; 
        0 1 1 0 1 1 1 0 0 0 0 0 1 0 0; 
        1 1 1 0 0 1 1 0 0 0 0 0 0 1 0; 
        1 0 1 0 0 0 1 0 0 0 0 0 0 0 1
    ];

    % 生成矩阵
    handles.G = [
        1 0 0 0 0 0 0 1 0 0 0 1 0 1 1;
        0 1 0 0 0 0 0 1 1 0 0 1 1 1 0;
        0 0 1 0 0 0 0 0 1 1 0 0 1 1 1;
        0 0 0 1 0 0 0 1 0 1 1 1 0 0 0;
        0 0 0 0 1 0 0 0 1 0 1 1 1 0 0;
        0 0 0 0 0 1 0 0 0 1 0 1 1 1 0;
        0 0 0 0 0 0 1 0 0 0 1 0 1 1 1
    ];

    % Update handles structure
    guidata(hObject, handles);
% --- Outputs from this function are returned to the command line.
function varargout = untitled4_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 打开文件选择对话框，选择图像文件
    [file, path] = uigetfile({'*.png;*.jpg;*.jpeg;*.bmp', 'Image Files (*.png, *.jpg, *.jpeg, *.bmp)'}, 'Select an Image');

    % 如果用户选择了图像文件
    if isequal(file, 0)
        disp('User canceled image selection.');
    else
        % 读取图像
        img = imread(fullfile(path, file));
        
        % 获取图像显示区域的句柄
        ax = handles.axes1; % 假设你在 GUI 中创建了一个名为 axes1 的图像显示区域
        
        % 在图像显示区域显示图像
        imshow(img, 'Parent', ax);
        title(ax, 'Loaded Image');
        
     handles.img = img;  % 将图像存入 handles 结构
     
     guidata(hObject, handles);  % 更新 handles
    end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    encoded =  handles.encoded ;
    nframe_num = handles.nframe_num;
    L = handles.L;
    Fs = handles.Fs;

    nframe = 0; % 当前帧
    jidai = [];
    while nframe < nframe_num
        nframe = nframe + 1;
        % 产生单极性波形
        fz = ones(1, L);
        msg = encoded(nframe, :) ;
        x1 = msg(fz, :);
        jidai(nframe, :) = reshape(x1, 1, L * length(msg)); % 基带信号
    
    end

    handles.jidai = jidai;  % 将图像存入 handles 结构
    guidata(hObject, handles);  % 更新 handles

    % 绘制时域图
    axes(handles.axes4);  % 指定时域图显示区域
    plot(jidai);  % 绘制时域图
    title('Waveform (Time Domain)');  
    xlabel('Sample Index');
    ylabel('Amplitude');

    % 绘制频谱图
    axes(handles.axes5);  % 指定频谱图显示区域
    % 计算频谱
    N = length(jidai);
    Y = fft(jidai);
    f = (0:N-1)*(Fs/N);  % 频率向量
    plot(f, abs(Y));  % 绘制频谱图
    title('Spectrum (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    jidai =  handles.jidai ;
    nframe_num = handles.nframe_num;
    Fs = handles.Fs;
    t = handles.t;
    nframe = 0; % 当前帧
    ask2 = [];
    while nframe < nframe_num
        nframe = nframe + 1;
        fc = 10000; % 载波频率
        zb = cos(2*pi*fc*t); % 载波
        msg = jidai(nframe, :) ;
        ask2(nframe, :)  = msg .* zb; % 2ASK调制
    end

    handles.ask2 = ask2;  % 将图像存入 handles 结构
    handles.zb = zb; 
    guidata(hObject, handles);  % 更新 handles

    % 绘制时域图
    axes(handles.axes4);  % 指定时域图显示区域
    plot(ask2);  % 绘制时域图
    title('Waveform (Time Domain)');  
    xlabel('Sample Index');
    ylabel('Amplitude');

    % 绘制频谱图
    axes(handles.axes5);  % 指定频谱图显示区域
    % 计算频谱
    N = length(ask2);
    Y = fft(ask2);
    f = (0:N-1)*(Fs/N);  % 频率向量
    plot(f, abs(Y));  % 绘制频谱图
    title('Spectrum (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
choose = get(handles.popupmenu7,'Value');

    db=get(handles.edit1,'String');
    db = str2double(db);  % 将字符串转换为 double 类型

    ask2 =  handles.ask2 ;
    Fs = handles.Fs;
    
    signal_power = mean(abs(ask2).^2);

     % 根据目标 SNR 计算噪声功率
    snr_linear = 10^(db / 10);  % 将 SNR 从 dB 转换为线性比例
    noise_power = signal_power / snr_linear;  % 噪声功率

    % 生成白噪声
    noise = sqrt(noise_power) .* randn(size(ask2));  % 生成高斯白噪声
    tz = awgn(ask2,db);
   

    handles.tz = tz;  % 将图像存入 handles 结构
    guidata(hObject, handles);  % 更新 handles

    % 绘制时域图
    axes(handles.axes4);  % 指定时域图显示区域
    plot(tz);  % 绘制时域图
    title('Waveform (Time Domain)');  
    xlabel('Sample Index');
    ylabel('Amplitude');

if choose == 1
    % 绘制频谱图
    axes(handles.axes5);  % 指定频谱图显示区域
    % 计算频谱
    N = length(tz);
    Y = fft(tz);
    f = (0:N-1)*(Fs/N);  % 频率向量
    plot(f, abs(Y));  % 绘制频谱图
    title('Spectrum (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
 else 
    % 绘制功率谱密度图
    axes(handles.axes5);  % 指定功率谱密度图显示区域
    
    N = length(tz);
    % 计算信号的功率谱密度
    % 使用 pwelch 函数计算功率谱密度，选择窗函数（例如Hamming窗）和重叠部分的长度
    [pxx, f] = pwelch(tz, hamming(N), N/2, N, Fs);
    
    % 绘制功率谱密度图
    plot(f, 10*log10(pxx));  % 转换为对数尺度（dB）
    title('Power Spectral Density (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');

 end


% 绘制功率谱密度图
plot(f, 10*log10(pxx));  % 转换为对数尺度（dB）
title('Power Spectral Density (Frequency Domain)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');


%     % 绘制噪声信号时域图
%     axes(handles.axes2);  % 指定 axes2 显示区域
%     plot(noise);  % 绘制加噪后的信号时域图
%     title(['Noisy Signal at ' num2str(db) ' dB SNR']);
%     xlabel('Sample Index');
%     ylabel('Amplitude');
    
    % 更新用户输入的信噪比显示
    disp(['The signal-to-noise ratio (SNR) is set to ' num2str(db) ' dB.']);

   

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double



% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    choose = get(handles.popupmenu7,'Value');

    tz =  handles.tz ;
    zb = handles.zb;
    Fs = handles.Fs;
   

   

         % 带通滤波器
        fp1 = 9000; fp2 = 11000; 
        fs1 = 1000; fs2 = 30000;
        wp1 = 2*fp1/Fs; wp2 = 2*fp2/Fs;
        ws1 = 2*fs1/Fs; ws2 = 2*fs2/Fs;
        % 1：通带最大衰减（dB）。
        % 30：阻带最小衰减（dB）。

        [N, wc] = buttord([wp1, wp2], [ws1, ws2], 1, 30); % 巴特沃斯滤波器
        [b, a] = butter(N, wc, 'bandpass');
        tz = filter(b, a, tz);
    
        % 相干解调
        tz = tz .* zb * 2; % 解调
    
        % 低通滤波
        fp = 6000; fs = 10000;
        wp = 2*fp/Fs; ws = 2*fs/Fs;
        [N, wc] = buttord(wp, ws, 1, 30);
        [b, a] = butter(N, wc, 'low');
        lvbo = filter(b, a, tz);

    handles.lvbo = lvbo;  % 将图像存入 handles 结构
    guidata(hObject, handles);  % 更新 handles

    % 绘制时域图
    axes(handles.axes4);  % 指定时域图显示区域
    plot(lvbo);  % 绘制时域图
    title('Waveform (Time Domain)');  
    xlabel('Sample Index');
    ylabel('Amplitude');

if choose == 1
    % 绘制频谱图
    axes(handles.axes5);  % 指定频谱图显示区域
    % 计算频谱
    N = length(lvbo);
    Y = fft(lvbo);
    f = (0:N-1)*(Fs/N);  % 频率向量
    plot(f, abs(Y));  % 绘制频谱图
    title('Spectrum (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
 else 
    % 绘制功率谱密度图
    axes(handles.axes5);  % 指定功率谱密度图显示区域
    
    N = length(lvbo);
    % 计算信号的功率谱密度
    % 使用 pwelch 函数计算功率谱密度，选择窗函数（例如Hamming窗）和重叠部分的长度
    [pxx, f] = pwelch(lvbo, hamming(N), N/2, N, Fs);
    
    % 绘制功率谱密度图
    plot(f, 10*log10(pxx));  % 转换为对数尺度（dB）
    title('Power Spectral Density (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');

 end

%     % 绘制带通滤波器频率响应
%     axes(handles.axes4);  % 在 axes1 上绘制
%     freqz(b, a, 1024, Fs);  % 绘制频率响应
%     title('Bandpass Filter Frequency Response');
%     xlabel('Frequency (Hz)');
%     ylabel('Magnitude (dB)');
% 
%     % 绘制低通滤波器频率响应
%     axes(handles.axes5);  % 在 axes2 上绘制
%     freqz(b, a, 1024, Fs);  % 绘制频率响应
%     title('Lowpass Filter Frequency Response');
%     xlabel('Frequency (Hz)');
%     ylabel('Magnitude (dB)');
% 绘制带通滤波器频率响应





function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    lvbo =  handles.lvbo ;
    nframe_num = handles.nframe_num;
    Fs = handles.Fs;

    nframe = 0; % 当前帧
    rx_code = [];
    while nframe < nframe_num
        nframe = nframe + 1;

        % 抽样判决
        rx_data = reshape(lvbo(nframe, :), 100, []); % 重构为 100x15 矩阵
%         rx_code(nframe, :) = mean(rx_data) > 0.5; % 判决
        rx_code(nframe, :) = rx_data(1, :) > 0.5; % 使用第一行数据进行判决
    end

    

    handles.rx_code = rx_code;  % 将图像存入 handles 结构
    guidata(hObject, handles);  % 更新 handles

      % 绘制时域图
    axes(handles.axes4);  % 指定时域图显示区域
    plot(rx_code);  % 绘制时域图
    title('Waveform (Time Domain)');  
    xlabel('Sample Index');
    ylabel('Amplitude');

    % 绘制频谱图
    axes(handles.axes5);  % 指定频谱图显示区域
    % 计算频谱
    N = length(rx_code);
    Y = fft(rx_code);
    f = (0:N-1)*(Fs/N);  % 频率向量
    plot(f, abs(Y));  % 绘制频谱图
    title('Spectrum (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
choose = get(handles.popupmenu7,'Value');

    rx_code_sequence =  handles.rx_code_sequence ;
    wave = handles.wave;
    Fs = handles.Fs;
    fz = handles.fz;
    n = handles.n;
    k = handles.k;
    H = handles.H;
    L = handles.L;

    
    

    decoded_sequence = Block_decoder(n, k, rx_code_sequence, H, 5);
    decoded_sequence =  decoded_sequence(:,1:7);
    fprintf('信道解码后的信号为: %s\n', num2str(decoded_sequence));


    x1=decoded_sequence(fz,:); 
    decoded=reshape(x1,1,L*k);    % 产生单极性不归零矩形脉冲波形，将刚得到的L*M矩阵，按列重新排列形成1*(L*M)的矩阵
    
    err = sum(wave ~= decoded_sequence);
    

    % 计算总的误码率
    total_errors = err;  % 误码总数
    total_bits = k;  % 总比特数（假设每帧有 n 个比特）
    total_ber = total_errors / total_bits;  % 总误码率
    
    % 更新GUI中显示总误码率
    set(handles.text19, 'String', ['Total BER: ', num2str(total_ber)]);
    fprintf('误码率 = %f \n', total_ber);


% %     receive_code = decoded(1,:);
%     new_row = mean(decoded,1);
%     receive_code = new_row';
%     disp(receive_code)
%     receive_code = receive_code(1,7);

% 
%     err= sum(receive_code ~= decoded);
%     fprintf('误码率 = %f \n', err);

%     nframe = 0; % 当前帧
%     receive_data = [];
%     while nframe < nframe_num
%         nframe = nframe + 1;
%         decoded_msg = Block_decoder(n, k, rx_code(nframe, :), H, 5);
%         receive_data(nframe, 1:7) = decoded_msg(1, 1:7);
% 
%          % 计算误码率
%         err(nframe) = sum(decoded_msg ~= encoded(nframe, :));
%         fprintf('误码率 = %f \n', err);
%     
%         if err
%             nferr = nferr + 1; % 误帧数
%         end
%     end
    handles.decoded = decoded;  % 将图像存入 handles 结构
    guidata(hObject, handles);  % 更新 handles
    
%     % 计算总的误码率
%     total_errors = sum(err);  % 误码总数
%     total_bits = nframe_num * 15;  % 总比特数（假设每帧有 n 个比特）
%     total_ber = total_errors / total_bits;  % 总误码率
    
%     % 更新GUI中显示总误码率
%     set(handles.text19, 'String', ['Total BER: ', num2str(err)]);

      % 绘制时域图
    axes(handles.axes4);  % 指定时域图显示区域
    plot(decoded);  % 绘制时域图
    title('Waveform (Time Domain)');  
    xlabel('Sample Index');
    ylabel('Amplitude');

   if choose == 1
    % 绘制频谱图
    axes(handles.axes5);  % 指定频谱图显示区域
    % 计算频谱
    N = length(decoded);
    Y = fft(decoded);
    f = (0:N-1)*(Fs/N);  % 频率向量
    plot(f, abs(Y));  % 绘制频谱图
    title('Spectrum (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
 else 
    % 绘制功率谱密度图
    axes(handles.axes5);  % 指定功率谱密度图显示区域
    
    N = length(decoded);
    % 计算信号的功率谱密度
    % 使用 pwelch 函数计算功率谱密度，选择窗函数（例如Hamming窗）和重叠部分的长度
    [pxx, f] = pwelch(decoded, hamming(N), N/2, N, Fs);
    
    % 绘制功率谱密度图
    plot(f, 10*log10(pxx));  % 转换为对数尺度（dB）
    title('Power Spectral Density (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');

 end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

receive_data =  handles.receive_data;
data_length = handles.data_length;

% 恢复数据流
rx_data_stream = reshape(receive_data', 1, []);
rx_data_stream = rx_data_stream(1:data_length);

% 恢复图像
rx_data_bin = num2str(rx_data_stream);
rx_data_bin = strrep(rx_data_bin, ' ', '');  % 去掉空格
img_bin = reshape(rx_data_bin, 8, [])'; 
img_reconstructed = bin2dec(img_bin);  
img_reconstructed = reshape(img_reconstructed, 32, 32); 

axes(handles.axes2);  % 指定axes2为显示区域
imshow(uint8(img_reconstructed));  % 显示恢复的图像
title('Reconstructed Image');


function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
 choose = get(handles.popupmenu7,'Value');

 
 n = handles.n;
 k = handles.k;
 G = handles.G;
 Fs = handles.Fs;
 L = handles.L;
  
 wave = handles.wave;
 fz = handles.fz;

 encoded_sequence = Block_encoder(n, k, wave, G); 
 fprintf('信道编码后的信号为: %s\n', num2str(encoded_sequence));

 encoded = encoded_sequence(fz,:);
 encoded = reshape(encoded,1,L*n);
 
 handles.encoded = encoded;
 handles.encoded_sequence = encoded_sequence;
 guidata(hObject, handles);  % 更新 handles

% 绘制时域图
    axes(handles.axes4);  % 指定时域图显示区域
    plot(encoded);  % 绘制时域图
    title('Waveform (Time Domain)');  
    xlabel('Sample Index');
    ylabel('Amplitude');

    if choose == 1
    % 绘制频谱图
    axes(handles.axes5);  % 指定频谱图显示区域
    % 计算频谱
    N = length(encoded);
    Y = fft(encoded);
    f = (0:N-1)*(Fs/N);  % 频率向量
    plot(f, abs(Y));  % 绘制频谱图
    title('Spectrum (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
 else 
    % 绘制功率谱密度图
    axes(handles.axes5);  % 指定功率谱密度图显示区域
    
    N = length(encoded);
    % 计算信号的功率谱密度
    % 使用 pwelch 函数计算功率谱密度，选择窗函数（例如Hamming窗）和重叠部分的长度
    [pxx, f] = pwelch(encoded, hamming(N), N/2, N, Fs);
    
    % 绘制功率谱密度图
    plot(f, 10*log10(pxx));  % 转换为对数尺度（dB）
    title('Power Spectral Density (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');

 end



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
    tiaozhi = get(handles.popupmenu4,'Value');
    choose = get(handles.popupmenu7,'Value');


    encoded =  handles.encoded ;
    Fs = handles.Fs;
    t = handles.t;
   
        switch tiaozhi
        case 1                        % 选中第一行
          fc = 10000; % 载波频率
          zb = cos(2*pi*fc*t); % 载波
          ask2 = encoded.* zb; % 2ASK调制
        case 2                        % 选中第二行
          
        case 3                        % 选中第三行
         
     
        end
        
    

    handles.ask2 = ask2;  % 将图像存入 handles 结构
    handles.zb = zb; 
    handles.tiaozhi = tiaozhi; 

    guidata(hObject, handles);  % 更新 handles

    % 绘制时域图
    axes(handles.axes4);  % 指定时域图显示区域
    plot(ask2);  % 绘制时域图
    title('Waveform (Time Domain)');  
    xlabel('Sample Index');
    ylabel('Amplitude');

    

 if choose == 1
    % 绘制频谱图
    axes(handles.axes5);  % 指定频谱图显示区域
    % 计算频谱
    N = length(ask2);
    Y = fft(ask2);
    f = (0:N-1)*(Fs/N);  % 频率向量
    plot(f, abs(Y));  % 绘制频谱图
    title('Spectrum (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
 else 
    % 绘制功率谱密度图
    axes(handles.axes5);  % 指定功率谱密度图显示区域
    
    N = length(ask2);
    % 计算信号的功率谱密度
    % 使用 pwelch 函数计算功率谱密度，选择窗函数（例如Hamming窗）和重叠部分的长度
    [pxx, f] = pwelch(ask2, hamming(N), N/2, N, Fs);
    
    % 绘制功率谱密度图
    plot(f, 10*log10(pxx));  % 转换为对数尺度（dB）
    title('Power Spectral Density (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');

 end

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    lvbo =  handles.lvbo ;
%     disp(lvbo);
    Fs = handles.Fs;
    fz = handles.fz;
    L = handles.L;
    n = handles.n;
    % 定义取样间隔
    interval = 100;
    % 获取取样点索引
    indices = 50:interval:length(lvbo);
    
    % 对取样点进行判断并生成结果列表
    rx_code_sequence = 1 * (lvbo(indices) > 0.5);

    x1=rx_code_sequence(fz,:); 
    rx_code=reshape(x1,1,L*n);    % 产生单极性不归零矩形脉冲波形，将刚得到的L*M矩阵，按列重新排列形成1*(L*M)的矩阵
    
    fprintf('抽样判决后的信号为: %s\n', num2str(rx_code_sequence));


    handles.rx_code_sequence = rx_code_sequence; 
    guidata(hObject, handles);  % 更新 handles

      % 绘制时域图
    axes(handles.axes4);  % 指定时域图显示区域
    plot(rx_code);  % 绘制时域图
    title('Waveform (Time Domain)');  
    xlabel('Sample Index');
    ylabel('Amplitude');

    % 绘制频谱图
    axes(handles.axes5);  % 指定频谱图显示区域
    % 计算频谱
    N = length(rx_code);
    Y = fft(rx_code);
    f = (0:N-1)*(Fs/N);  % 频率向量
    plot(f, abs(Y));  % 绘制频谱图
    title('Spectrum (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

choose = get(handles.popupmenu7,'Value');

    k = handles.k;
    Fs = handles.Fs;
    L = handles.L;

    
%     wave=randi([0,1],1,k); % 产生的随机信号
    wave = [0,1,0,0,1,1,1];
    fz=ones(1,L);
    x1=wave(fz,:); 
    jidai=reshape(x1,1,L*k);    % 产生单极性不归零矩形脉冲波形，将刚得到的L*M矩阵，按列重新排列形成1*(L*M)的矩阵
    
    handles.wave = wave;
    handles.x1 = x1;
    handles.fz = fz;
    handles.jidai = jidai;

    fprintf('产生的二进制随机序列为: %s\n', num2str(wave));
    guidata(hObject, handles);  % 更新 handles

    % 绘制时域图
    axes(handles.axes4);  % 指定时域图显示区域
    plot(jidai);  % 绘制时域图
    title('Waveform (Time Domain)');  
    xlabel('Sample Index');
    ylabel('Amplitude');

 if choose == 1
    % 绘制频谱图
    axes(handles.axes5);  % 指定频谱图显示区域
    % 计算频谱
    N = length(jidai);
    Y = fft(jidai);
    f = (0:N-1)*(Fs/N);  % 频率向量
    plot(f, abs(Y));  % 绘制频谱图
    title('Spectrum (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
 else 
    % 绘制功率谱密度图
    axes(handles.axes5);  % 指定功率谱密度图显示区域
    
    N = length(jidai);
    % 计算信号的功率谱密度
    % 使用 pwelch 函数计算功率谱密度，选择窗函数（例如Hamming窗）和重叠部分的长度
    [pxx, f] = pwelch(jidai, hamming(N), N/2, N, Fs);
    
    % 绘制功率谱密度图
    plot(f, 10*log10(pxx));  % 转换为对数尺度（dB）
    title('Power Spectral Density (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');

 end

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    db=get(handles.edit1,'String');
    db = str2double(db);  % 将字符串转换为 double 类型

    img = handles.img;
    n = handles.n;
    k = handles.k;
    L = handles.L;
    G = handles.G;
    H = handles.H;
    Fs = handles.Fs;
    t = handles.t;

    w=120;h=120;
    input_image=reshape(img,1,w*h,3);
    input_data = reshape(sum(permute(input_image, [3 1 2])) > 255*3/2, size(input_image, 1), size(input_image, 2));

    data_length = length(input_data);
    pad_length = mod(-data_length, k);
    
    input_data_padded = [input_data, zeros(1, pad_length)];
    grouped_data = reshape(input_data_padded, k, [])';
    nframe_num = size(grouped_data, 1); %帧数
    
    encoded_msg = zeros(nframe_num, n);
    
    nframe = 1; % 当前帧
    fz = ones(1, L);
    fc = 10000; % 载波频率
    zb = cos(2*pi*fc*t); % 载波
    receive_data = [];
    err = [];
    errs = 0;
    nferr = 0;
    while nframe <= nframe_num
        
        msg = grouped_data(nframe, :); % 获取每一行信息
        encoded_msg(nframe, :) = Block_encoder(n, k, msg, G); % 编码
        
        msg = encoded_msg(nframe, :);
        x1 = msg(fz, :);
        msg = reshape(x1, 1, L * length(msg)); % 基带信号
        msg  = msg .* zb; % 2ASK调制
        msg  = awgn(msg,db);


         % 带通滤波器
        fp1 = 9000; fp2 = 11000; 
        fs1 = 1000; fs2 = 30000;
        wp1 = 2*fp1/Fs; wp2 = 2*fp2/Fs;
        ws1 = 2*fs1/Fs; ws2 = 2*fs2/Fs;
        [N, wc] = buttord([wp1, wp2], [ws1, ws2], 1, 30);
        [b, a] = butter(N, wc, 'bandpass');
        msg = filter(b, a, msg);
    
        % 相干解调
        msg = msg .* zb * 2; % 解调
    
        % 低通滤波
        fp = 6000; fs = 10000;
        wp = 2*fp/Fs; ws = 2*fs/Fs;
        [N, wc] = buttord(wp, ws, 1, 30);
        [b, a] = butter(N, wc, 'low');
        msg = filter(b, a, msg);

        msg = reshape(msg, 100, []); % 重构为 100x15 矩阵
        msg = msg(50,:) > 0.5;
%         msg = mean(msg) > 0.5; % 判决
        decoded_msg = Block_decoder(n, k, msg, H, 5);
        receive_data(nframe, 1:7) = decoded_msg(1, 1:7);

        err(nframe) = sum(decoded_msg ~= encoded_msg(nframe, :));
 
        fprintf('误码率 = %f \n', err(nframe)/15);
        
    
        if err(nframe)
            nferr = nferr + 1; % 误帧数
        end
        nframe = nframe + 1;
    
    end
        
        disp("图像传输成功")
        fprintf('总误码数 = %f \n', sum(err));
        fprintf('总误码率 = %f \n', mean(err)/15);
        fprintf('误帧数 = %f \n', nferr);
        fprintf('误帧率 = %f \n', nferr/nframe_num);

    % 图片解压
%     output = repmat(receive_data * 255, 1, 3); % 数据恢复到像素值范围
        receive_data_padded = reshape(receive_data', 1, []); % 转置并展开
        data_length = w * h; % 原始数据长度
        receive_data_recovered = receive_data_padded(1:data_length); % 移除填充

        % 2. 将逻辑矩阵还原为原始图片的二值化形式
        receive_image_recovered = reshape(receive_data_recovered, h, w); % 还原到二维图片

        % 3. 重构原始 RGB 图像
        img_recovered = uint8(repmat(receive_image_recovered, 1, 1, 3) * 255); % 二值化转为 RGB 格式

     % 获取图像显示区域的句柄
        ax = handles.axes2; % 假设你在 GUI 中创建了一个名为 axes1 的图像显示区域
        
        % 在图像显示区域显示图像
        imshow(img_recovered, 'Parent', ax);
        title(ax, 'Loaded Image');
% %     绘制误码率图
% %     axes(handles.axes2);  % 在axes2上绘制误码率
% %     plot(err);  % 绘制误码数图像
% %     title('Bit Error Rate (BER) per Frame');
% %     xlabel('Frame Index');
% %     ylabel('Bit Error Rate');



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double

choose = get(handles.popupmenu7,'Value');

    db=get(handles.edit9,'String');
    threshould = str2double(db);  % 将字符串转换为 double 类型

    lvbo =  handles.lvbo ;
%     disp(lvbo);
    Fs = handles.Fs;
    fz = handles.fz;
    L = handles.L;
    n = handles.n;
    % 定义取样间隔
    interval = 100;
    % 获取取样点索引
    indices = 50:interval:length(lvbo);
    
    % 对取样点进行判断并生成结果列表
    rx_code_sequence = 1 * (lvbo(indices) > threshould);

    x1=rx_code_sequence(fz,:); 
    rx_code=reshape(x1,1,L*n);    % 产生单极性不归零矩形脉冲波形，将刚得到的L*M矩阵，按列重新排列形成1*(L*M)的矩阵
    
    fprintf('抽样判决后的信号为: %s\n', num2str(rx_code_sequence));


    handles.rx_code_sequence = rx_code_sequence; 
    guidata(hObject, handles);  % 更新 handles

      % 绘制时域图
    axes(handles.axes4);  % 指定时域图显示区域
    plot(rx_code);  % 绘制时域图
    title('Waveform (Time Domain)');  
    xlabel('Sample Index');
    ylabel('Amplitude');

if choose == 1
    % 绘制频谱图
    axes(handles.axes5);  % 指定频谱图显示区域
    % 计算频谱
    N = length(rx_code);
    Y = fft(rx_code);
    f = (0:N-1)*(Fs/N);  % 频率向量
    plot(f, abs(Y));  % 绘制频谱图
    title('Spectrum (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
 else 
    % 绘制功率谱密度图
    axes(handles.axes5);  % 指定功率谱密度图显示区域
    
    N = length(rx_code);
    % 计算信号的功率谱密度
    % 使用 pwelch 函数计算功率谱密度，选择窗函数（例如Hamming窗）和重叠部分的长度
    [pxx, f] = pwelch(rx_code, hamming(N), N/2, N, Fs);
    
    % 绘制功率谱密度图
    plot(f, 10*log10(pxx));  % 转换为对数尺度（dB）
    title('Power Spectral Density (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');

 end


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbutton13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
