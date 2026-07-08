%% 第一步：配置必要的前置参数（Trans/Receive/Resource 结构体）
% 1. 发射参数结构体（Trans）：定义发射频率
Trans.frequency = 7.6;  % 发射频率 7.6 MHz（L11-5  transducer 常用频率）

% 2. 接收参数结构体（Receive）：定义采样模式、抽取率等核心参数
% 配置 2 套滤波器（演示多配置集功能）
Receive(1).sampleMode = 'NS200BW';        % 采样模式：4样本/波长，200%带宽
Receive(1).decimSampleRate = 4 * Trans.frequency;  % 抽取采样率 = 4×发射频率（默认匹配）
Receive(1).demodFrequency = Trans.frequency;       % 解调频率 = 发射频率
Receive(1).LowPassCoef = [];  % 初始为空，将由 filterTool 生成
Receive(1).InputFilter = [];  % 初始为空，将由 filterTool 生成

Receive(2).sampleMode = 'BS100BW';        % 第二套配置：2样本/波长，100%带宽
Receive(2).decimSampleRate = 2 * Trans.frequency;  % 抽取采样率 = 2×发射频率
Receive(2).demodFrequency = Trans.frequency;
Receive(2).LowPassCoef = [];
Receive(2).InputFilter = [];

% 3. 系统资源结构体（Resource）：设置仿真模式（可选，避免实时模式限制）
Resource.Parameters.simulateMode = 0;  % 0=实时模式，1=仿真模式（仅可视化不联动硬件）

% 将结构体存入工作区（filterTool 需读取这些变量）
assignin('base', 'Trans', Trans);
assignin('base', 'Receive', Receive);
assignin('base', 'Resource', Resource);

%% 第二步：启动 filterTool 滤波器设计工具
filterTool;  % 无输入参数 = 常规模式，自动读取工作区变量

%% 第三步：（可选）设计完成后，提取滤波器系数并应用到系统
% 工具启动后，手动调整参数并点击「Filter Generation」，命令行会输出系数
% 以下是示例：将生成的系数更新到 Receive 结构体，用于后续成像
% Receive(1).LowPassCoef = [+0.0000 +0.0000 +0.0000 +0.0000 +0.0000 +0.0000 ...
%                          +0.0000 +0.0000 +0.0000 +0.0000 +0.0000 +1.0000];  % 示例系数
% Receive(1).InputFilter = [-0.0011 +0.0000 -0.0012 +0.0000 +0.0055 +0.0000 ...
%                          +0.0072 +0.0000 -0.0142 +0.0000 -0.0264 +0.0000 ...
%                          +0.0261 +0.0000 +0.0782 +0.0000 -0.0367 +0.0000 ...
%                          -0.3079 +0.0000 +0.5411];  % 示例系数

% % 更新工作区的 Receive 结构体
% assignin('base', 'Receive', Receive);
% 
% % （系统环境）联动 Vantage 系统更新滤波器
% if evalin('base', 'exist(''Control'',''var'')')
%     Control = evalin('base', 'Control');
%     Control(length(Control)+1).Command = 'update&Run';
%     Control(end).Parameters = {'Receive'};
%     assignin('base', 'Control', Control);
% end