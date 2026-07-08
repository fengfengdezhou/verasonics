% =========================================================================
% 名称：Stimulate_Once_Then_Image_GUI.m
% 功能：(V2) 执行一次GUI控制的聚焦刺激, 然后进入连续的成像循环以观察变化。
%
% 序列流程:
%   - GUI 提供两个启动按钮:
%     1. "开始成像" (startEvent = 3):
%        - 系统进入 [成像 -> 重建 -> 跳转回成像] 的循环 (Event 3 -> N -> 3)
%     2. "刺激并成像" (startEvent = 1):
%        - 系统执行 Event 1 (刺激)
%        - Event 2 (跳转到 Event 3)
%        - 系统进入 [成像 -> 重建 -> 跳转回成像] 的循环 (Event 3 -> N -> 3)
%
% GUI功能:
%   - 滑块实时调节刺激脉冲 (TX(na+1)) 的焦点 (X 和 Z)。
%   - 按钮可触发 "刺激并成像" 或 "仅成像" 序列。
%   - 按钮可单独预览声场而不运行序列。
% =========================================================================

clear all; clc; close all;

% --- 基本参数定义 ---
P.startDepth = 15; % 成像开始深度 (wavelengths)
P.endDepth = 135;  % 成像结束深度 (wavelengths)
numRcvElements = 64; % 定义接收通道数（适配64LE）

% --- 多角度复合参数 ---
na = 3; % 角度数量
dtheta = (12 * pi / 180) / (na - 1);
P.startAngle = -12 * pi / 180 / 2;

% --- 1. 定义系统参数 (Resource) ---
Resource.VDAS.dmaTimeout = 10000;
Resource.Parameters.numTransmit = 128;   % 发射通道数
Resource.Parameters.numRcvChannels = 64;   % 接收通道数
Resource.Parameters.speedOfSound = 1540;   % 声速 (m/sec)
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 1; % 1 = 仿真模式, 0 = 运行硬件

% <<< 关键：设置默认启动事件为 "仅成像" 循环的入口 (Event 3) >>>
Resource.Parameters.startEvent = 3; 

% --- 2. 定义换能器 (Trans) ---
Trans.name = 'L14-6';
Trans.units = 'mm';
Trans.numelements = 128;
Trans.frequency = 9.25; % 探头中心频率 (MHz)
Trans.Bandwidth = [9.25*(1-1/2), 9.25*(1+1/2)];
Trans.type = 0; % 线性探头
Trans.id = hex2dec('0000');
Trans.connType = 1; % 接头类型
Trans.elevationApertureMm = 4;
Trans.elevationFocusMm = 1000000;
Trans.elementWidth = 0.16; 
Trans.spacingMm = 0.2; 
Trans.ElementPos = zeros(Trans.numelements,4);
Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
if ~isfield(Trans,'ElementSens')
    Theta = (-pi/2:pi/100:pi/2);
    Theta(51) = 0.0000001; 
    eleWidthWl = Trans.elementWidth * Trans.frequency/Resource.Parameters.speedOfSound;
    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
end
Trans.lensCorrection = 0;
Trans.impedance = 50;
Trans.maxHighVoltage = 30; 
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000); 
Trans.spacing = Trans.spacingMm * scaleToWvl; 

% --- 3. 定义像素数据 (PData) ---
PData(1).PDelta = [1, 0, 1] * 0.5;    
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); 
PData(1).Size(2) = ceil((numRcvElements * Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1; 
PData(1).Origin = [-Trans.spacing*(numRcvElements - 1)/2,0,P.startDepth];

% --- 4. 定义仿真介质 (Media) ---
Media.MP(1,:) = [0,0,50,1.0];
Media.MP(2,:) = [-10,0,60,1.0];
Media.MP(3,:) = [10,0,70,1.0];
Media.function = 'movePoints'; % <<< 添加一个移动函数用于演示变化 >>>

% --- 5. 定义资源缓存 (Resource Buffers) ---
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = na*4096*2; 
Resource.RcvBuffer(1).colsPerFrame = numRcvElements;
Resource.RcvBuffer(1).numFrames = 30; 
Resource.InterBuffer(1).datatype = 'complex single';
Resource.InterBuffer(1).numFrames = 1; 
Resource.ImageBuffer(1).numFrames = 10; 
Resource.DisplayWindow(1).Title = 'L14-6 Stimulate-Once-Then-Image';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% --- 6. 定义发射波形 (TW) ---
% TW(1) 用于成像 (0.5 周期)
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,0.5,1];
% TW(2) 用于刺激 (例如: 5 周期 = 10 个半周期)
TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,0.67,10,1]; 

% --- 7. 定义发射事件 (TX) ---
% TX(1:na) 用于成像
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', kaiser(Resource.Parameters.numTransmit,1)', ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, na);
for n = 1 : na 
    TX(n).Steer = [(P.startAngle + (n - 1) * dtheta),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end

% TX(na+1) 用于刺激 (Push)
stimTX.waveform = 2; 
stimTX.Origin = [0,0,0]; 
stimTX.Apod = ones(1,Trans.numelements); 
stimTX.focus = 30 * scaleToWvl; % 默认 Z=30mm
stimTX.Steer = [0,0];
stimTX.Delay = computeTXDelays(stimTX);
TX(na+1) = stimTX;

% --- 8. 定义时间增益补偿 (TGC) ---
TGC.CntrlPts = [489,519,565,611,641,687,733,779];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% --- 9. 定义接收事件 (Receive) ---
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((numRcvElements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, na * Resource.RcvBuffer(1).numFrames);

for i = 1 : Resource.RcvBuffer(1).numFrames
    % <<< 关键：在成像循环的第一个采集上调用 Media.function >>>
    % 这将使得 *只有* 在成像循环中目标点才会移动，
    % 在按下 "Stimulate" 按钮的瞬间目标点是静止的。
    Receive(na * (i - 1) + 1).callMediaFunc = 1; 
    
    for j = 1 : na
        Receive(na * (i - 1) + j).Apod(33 : 96) = 1.0; 
        Receive(na * (i - 1) + j).framenum = i;
        Receive(na * (i - 1) + j).acqNum = j;
    end
end

% --- 10. 定义重建 (Recon) ---
Recon = struct('senscutoff', 0.7, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ... 
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ... 
               'RINums', 1 : na);

% --- 11. 定义重建信息 (ReconInfo) ---
ReconInfo = repmat(struct('mode', 'accumIQ', ...
                          'txnum', 1, ...
                          'rcvnum', 1, ...
                          'regionnum', 1), 1, na);
if na > 1
    ReconInfo(1).mode = 'replaceIQ'; 
    for j = 1 : na 
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
    end
    ReconInfo(na).mode = 'accumIQ_replaceIntensity'; 
else
    ReconInfo(1).mode = 'replaceIntensity'; 
end

% --- 12. 定义处理 (Process) ---
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...
                         'framenum',-1,...
                         'pdatanum',1,...
                         'pgain',1.0,...
                         'reject',2,...
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...
                         'displayWindow',1};

% --- 13. 定义序列控制 (SeqControl) ---
% SeqControl(1): 跳转到 "仅成像" 循环 (Event 3)
SeqControl(1).command = 'jump'; 
SeqControl(1).argument = 3; % <<< 修改：跳转到 Event 3 >>>

SeqControl(2).command = 'timeToNextAcq';  % 成像角度之间的延迟
SeqControl(2).argument = 160;  % 160 us

SeqControl(3).command = 'timeToNextAcq';  % 帧间延迟 (确保总帧率)
SeqControl(3).argument = 20000 - (na-1)*SeqControl(2).argument; % 20 ms (50 fps)

SeqControl(4).command = 'returnToMatlab'; % 返回 Matlab (用于GUI响应)

SeqControl(5).command = 'timeToNextAcq';  % 刺激(Stim) 和 成像(Image) 间的延迟
SeqControl(5).argument = 200;  % 200 us

% <<< 添加: SeqControl(6) 用于从刺激事件跳转到成像循环 >>>
SeqControl(6).command = 'jump';
SeqControl(6).argument = 3; % 跳转到 Event 3

nsc = 7; % SeqControl 计数器

% --- 14. 定义事件序列 (Event) ---
% (V2: 具有两个入口的事件列表)
n = 1;

% === 入口 1: 刺激 (Event 1-2) ===
% (仅当 Resource.Parameters.startEvent = 1 时执行)
Event(n).info = 'Stimulation Pulse';
Event(n).tx = na + 1; % 使用 TX(na+1) [刺激]
Event(n).rcv = 0; % 不接收
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 5; % SeqControl(5) = 200 us 刺激-成像延迟
n = n + 1;

Event(n).info = 'Jump to Imaging Loop';
Event(n).tx = 0; 
Event(n).rcv = 0; 
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 6; % SeqControl(6) = 跳转到 Event 3
n = n + 1;

% === 入口 2: 仅成像循环 (Event 3 开始) ===
% (默认 Resource.Parameters.startEvent = 3)
for i = 1 : Resource.RcvBuffer(1).numFrames
    
    % --- 成像采集 (na 个角度) ---
    for j = 1 : na  
        Event(n).info = 'Imaging Angle';
        Event(n).tx = j; % 使用 TX(j) [成像]
        Event(n).rcv = na * (i - 1) + j; % 对应的接收事件
        Event(n).recon = 0;
        Event(n).process = 0;
        if j < na
            Event(n).seqControl = 2; % 角度间延迟
        else
            % 最后一个角度，设置帧率并触发传输
            Event(n).seqControl = [3, nsc]; % 帧延迟
            SeqControl(nsc).command = 'transferToHost'; % 传输数据
            nsc = nsc + 1;
        end
        n = n + 1;
    end
    
    % --- 重建 & 处理 ---
    Event(n).info = 'Recon and Process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 4; % returnToMatlab
    n = n + 1;
end

% --- 序列结尾跳转 ---
Event(n).info = 'Jump back to Imaging Loop Start';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1; % SeqControl(1) = 跳转到 Event 3
n = n + 1;


% --- 15. 定义 GUI 控件 (UI) ---
% (成像控制)
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

MinMaxVal = [64,300,P.endDepth]; 
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% % (刺激控制)
% UI(3).Control = {'UserB3','Style','VsSlider','Label','Stim X Offset (mm)',...
%                  'SliderMinMaxVal',[-5,5,0],'SliderStep',[0.05,0.2],'ValueFormat','%2.1f','Tag','xSliderTag'};
% UI(3).Callback = text2cell('%StimFocusUpdateCallback'); % <<< 修改: 仅更新参数 >>>
% 
% UI(4).Control = {'UserB2','Style','VsSlider','Label','Stim Z Depth (mm)',...
%                  'SliderMinMaxVal',[10,60,30],'SliderStep',[0.1,0.2],'ValueFormat','%2.0f','Tag','zSliderTag'};
% UI(4).Callback = text2cell('%StimFocusUpdateCallback'); % <<< 修改: 仅更新参数 >>>

% (序列触发按钮)
UI(5).Control = {'UserC2','Style','VsPushButton','Label','Preview Stim Field'};
UI(5).Callback = text2cell('%ShowSingleCallback'); % (预览声场)

% <<< 修改: 按钮 6 用于触发 "刺激并成像" (startEvent=1) >>>
UI(6).Control = {'UserC1','Style','VsPushButton','Label','Fire STIM & Image'};
UI(6).Callback = text2cell('%FireStimCallback');

% <<< 添加: 按钮 7 用于触发 "仅成像" (startEvent=3) >>>
UI(7).Control = {'UserC3','Style','VsPushButton','Label','Start IMAGING Only'};
UI(7).Callback = text2cell('%StartImagingCallback');

% --- 16. 保存 .mat 文件并运行 VSX ---
save('MatFiles/Stimulate_Once_Then_Image_GUI');
filename = 'Stimulate_Once_Then_Image_GUI';
VSX;
return

% =========================================================================
% **** 回调函数 (Callback Routines)                   ****
% =========================================================================

%SensCutoffCallback - (来自 L14-6_multiangle_LE)
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutoffCallback

%RangeChangeCallback - (来自 L14-6_multiangle_LE)
simMode = evalin('base','Resource.Parameters.simulateMode');
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        P.endDepth = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);
evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((numRcvElements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback

% %StimFocusUpdateCallback
% % (V2) 仅更新滑块参数，不触发序列
% Control = evalin('base','Control');
% if strcmp(Control.Command, 'update&Run')
%     return
% end
% xOffset = get(findobj('Tag','xSliderTag'),'Value');
% zDepth  = get(findobj('Tag','zSliderTag'),'Value');
% na = evalin('base', 'na');
% Trans = evalin('base', 'Trans');
% Resource = evalin('base', 'Resource');
% TX = evalin('base', 'TX');
% scaleToWvl = Trans.frequency / (Resource.Parameters.speedOfSound/1000);
% TX(na+1).Origin = [xOffset, 0, 0];
% TX(na+1).focus = zDepth * scaleToWvl;
% TX(na+1).Delay = computeTXDelays(TX(na+1));
% assignin('base', 'TX', TX);
% Control.Command = 'update&Run';
% Control.Parameters = {'TX'};
% assignin('base','Control', Control);
% disp(['[Callback] 刺激焦点参数已更新: X = ',num2str(xOffset),' mm, Z = ',num2str(zDepth),' mm']);
% return
% %StimFocusUpdateCallback

%FireStimCallback
% (V2) "Fire STIM & Image" 按钮的回调
% 1. (可选) 确保焦点参数是最新
evalin('base', 'StimFocusUpdateCallback');
% 2. 设置 startEvent 为 1
disp('[Callback] 序列启动：将从 Event(1) [刺激] 开始');
Resource = evalin('base', 'Resource');
Resource.Parameters.startEvent = 1; % 刺激入口
assignin('base', 'Resource', Resource);
% 3. 触发 set&Run
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters', 1,'startEvent',1};
assignin('base','Control', Control);
return
%FireStimCallback

%StartImagingCallback
% (V2) "Start IMAGING Only" 按钮的回调
disp('[Callback] 序列启动：将从 Event(3) [仅成像] 开始');
Resource = evalin('base', 'Resource');
Resource.Parameters.startEvent = 3; % 仅成像入口
assignin('base', 'Resource', Resource);
% 2. 触发 set&Run
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1, 'startEvent',3};
assignin('base','Control', Control);
return
%StartImagingCallback

% （保留）声场预览回调函数（读取固定的TX(na+1)参数）
%ShowSingleCallback
% 'Preview Stim Field' (预览声场) 按钮的回调函数
disp('[Callback] 正在计算单点刺激声场预览...');
xOffset = 0; % 与固定的TX(na+1).Origin(1)一致（无需修改）
zDepth  = 30; % 与固定的TX(na+1).focus/scaleToWvl一致（无需修改）
na = evalin('base', 'na');
Trans = evalin('base', 'Trans');
Resource = evalin('base', 'Resource');
PData_img = evalin('base', 'PData(1)'); 
TX = evalin('base', 'TX');
PData_preview = PData_img;
scaleToWvl = Trans.frequency / (Resource.Parameters.speedOfSound/1000);
xHalfRange = 10; 
zHalfRange = 10; 
xRangeWvl = xHalfRange * scaleToWvl;
zStartWvl = max(0, (zDepth - zHalfRange) * scaleToWvl);
zEndWvl   = (zDepth + zHalfRange) * scaleToWvl;
PData_preview.Origin = [xOffset - xRangeWvl, 0, zStartWvl];
PData_preview.Size = [ceil((zEndWvl - zStartWvl)/PData_preview.PDelta(3)), ...
                      ceil((2*xRangeWvl)/PData_preview.PDelta(1)), 1];
TX_stim = TX(na+1); % 直接读取固定后的刺激TX参数
TX_stim.TXPD = computeTXPD(TX_stim, PData_preview);
assignin('base', 'TX_preview', TX_stim);
assignin('base', 'PData_preview', PData_preview);
% 正确代码：第二个参数传递 PData 结构体实体 PData_preview
showTXPD('TX_preview', PData_preview); % 生成预览图 % 生成预览图
disp('[Callback] 声场预览完成.');
return
%ShowSingleCallback

%movePoints
% (V2) 添加一个简单的仿真移动函数，以便观察刺激前后的变化
Media = evalin('base','Media');
Media.MP(1,3) = Media.MP(1,3) + 0.1; % 向下移动 0.1 mm
Media.MP(2,1) = Media.MP(2,1) + 0.05; % 向右移动 0.05 mm
assignin('base','Media',Media);
return
%movePoints