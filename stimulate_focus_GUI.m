function stimulate_focus_GUI()
clc; close all;

%% ---------------------- 系统与探头参数 ------------------------
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.numTransmit = 128;
Resource.Parameters.numRcvChannels = 64;
assignin('base','Resource',Resource);

Trans.name = 'L14-6';
Trans.units = 'mm';
Trans.numelements = 128;
Trans.frequency = 9.25;
Trans.Bandwidth = [Trans.frequency*(1-0.5), Trans.frequency*(1+0.5)];
Trans.type = 0;
Trans.elementWidth = 0.16;
Trans.spacingMm = 0.2;

Trans.ElementPos = zeros(Trans.numelements,4);
Trans.ElementPos(:,1) = Trans.spacingMm * ( -((Trans.numelements-1)/2) : ((Trans.numelements-1)/2) );
Trans.ElementPos(:,2:4) = 0;

Theta = (-pi/2:pi/100:pi/2);
Theta(51)=1e-7;
eleWidthWl = Trans.elementWidth*Trans.frequency/Resource.Parameters.speedOfSound;
Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
scaleToWvl = Trans.frequency / (Resource.Parameters.speedOfSound/1000);
Trans.spacing = Trans.spacingMm * scaleToWvl;
Trans.maxHighVoltage = 60;  % 探头最大电压，仿真可用 60V
assignin('base','Trans',Trans);

%% ---------------------- 声场分辨率 ------------------------
PData(1).PDelta = [0.25, 0, 0.25];  % 分辨率

%% ---------------------- 发射波形 ------------------------
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,1,1];
assignin('base','TW',TW);

TX(1).waveform = 1;
TX(1).Apod = ones(1,Trans.numelements);
TX(1).focus = 30;    
TX(1).Steer = [0,0];
TX(1).Delay = computeTXDelays(TX(1));
assignin('base','TX',TX);

%% ---------------------- GUI界面 ------------------------
f = figure('Name','多点聚焦刺激控制','Position',[200 200 480 300],'Color',[0.9 0.9 0.9]);

uicontrol('Style','text','Position',[40 250 100 20],'String','横向偏移 (mm)','BackgroundColor',[0.9 0.9 0.9]);
xSlider = uicontrol('Style','slider','Min',-5,'Max',5,'Value',0,'Position',[150 250 250 20],'Callback',@updateFocus);

uicontrol('Style','text','Position',[40 200 100 20],'String','聚焦深度 (mm)','BackgroundColor',[0.9 0.9 0.9]);
zSlider = uicontrol('Style','slider','Min',10,'Max',60,'Value',30,'Position',[150 200 250 20],'Callback',@updateFocus);

uicontrol('Style','pushbutton','Position',[100 130 120 40],'String','显示当前焦点声场','FontSize',10,'Callback',@showSingle);
uicontrol('Style','pushbutton','Position',[260 130 120 40],'String','执行多点扫描','FontSize',10,'Callback',@multiScan);

uicontrol('Style','text','Position',[100 70 300 20],'String','可通过滑条调节焦点，或自动扫描多个焦点','BackgroundColor',[0.9 0.9 0.9]);

%% ---------------------- 回调函数 ------------------------
    function updateFocus(~,~)
        xOffset = get(xSlider,'Value');
        zDepth  = get(zSlider,'Value');
        disp(['当前焦点位置: X = ',num2str(xOffset),' mm, Z = ',num2str(zDepth),' mm']);
    end

    function showSingle(~,~)
        xOffset = get(xSlider,'Value');
        zDepth  = get(zSlider,'Value');

        % 动态设置声场显示范围
        xHalfRange = 10;  % ±10 mm
        zHalfRange = 10;  

        scaleToWvl = Trans.frequency / (Resource.Parameters.speedOfSound/1000);
        xRangeWvl = xHalfRange * scaleToWvl;
        zStartWvl = max(0, (zDepth - zHalfRange) * scaleToWvl);
        zEndWvl   = (zDepth + zHalfRange) * scaleToWvl;

        PData(1).Origin = [-xRangeWvl,0,zStartWvl];
        PData(1).Size = [ceil((zEndWvl - zStartWvl)/PData(1).PDelta(3)), ...
                         ceil((2*xRangeWvl)/PData(1).PDelta(1)), 1];
        assignin('base','PData',PData);

        % 更新 TX
        TX(1).Origin = [xOffset,0,0];
        TX(1).focus = zDepth * scaleToWvl;
        TX(1).Delay = computeTXDelays(TX(1));
        assignin('base','TX',TX);

        disp('计算声场中...');
        TX(1).TXPD = computeTXPD(TX(1),PData(1));
        assignin('base','TX',TX);
        showTXPD(1);
    end

    function multiScan(~,~)
        disp('执行多点扫描...');
        xPoints = linspace(-3,3,5);
        zPoints = [25 35 45];
        for z = zPoints
            for x = xPoints
                % 动态设置声场显示范围
                xHalfRange = 10;  
                zHalfRange = 10;  
                scaleToWvl = Trans.frequency / (Resource.Parameters.speedOfSound/1000);
                xRangeWvl = xHalfRange * scaleToWvl;
                zStartWvl = max(0, (z - zHalfRange) * scaleToWvl);
                zEndWvl   = (z + zHalfRange) * scaleToWvl;

                PData(1).Origin = [-xRangeWvl,0,zStartWvl];
                PData(1).Size = [ceil((zEndWvl - zStartWvl)/PData(1).PDelta(3)), ...
                                 ceil((2*xRangeWvl)/PData(1).PDelta(1)), 1];
                assignin('base','PData',PData);

                TX(1).Origin = [x,0,0];
                TX(1).focus = z * scaleToWvl;
                TX(1).Delay = computeTXDelays(TX(1));
                TX(1).TXPD = computeTXPD(TX(1),PData(1));
                assignin('base','TX',TX);

                disp(['→ X=',num2str(x),' mm, Z=',num2str(z),' mm']);
                showTXPD(1);
                pause(0.5);
            end
        end
        disp('多点扫描完成！');
    end
end
