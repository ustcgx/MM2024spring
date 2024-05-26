% MATLAB仿真程序
clear;
clc;

%% 参数设置
Lx = 300; % x轴模拟范围 m
Ly = 300; % y轴模拟范围 m
Lz = 100; % z轴模拟范围 m
T = 1000000; % 模拟时间 s

figx = 200; % x轴可视化范围 m
figy = 200; % y轴可视化范围 m
figz = 100; % z轴可视化范围 m

dx = 1; % x方向网格大小 m
dy = 1; % y方向网格大小 m
dz = 1; % z方向网格大小 m
dt = 0.1; % 时间步长 s

D = 1; % 扩散系数 m^2/s，取相对大值
alpha = 1.0 % 反射系数，1.0代表全反射，0代表全吸收
v = 1; % 风速，沿x轴方向 m/s
w = 0.006; % 模拟烟尘重力的风速，沿z轴负向 m/s
C0 = 1e-2; % 烟囱排放的恒定浓度 kg/s
h = 20; % 烟囱高度 m

%% 创建网格
num_x = round(2*Lx/dx)+1;
num_y = round(2*Ly/dy)+1;
num_z = round(Lz/dz)+1;
[x, y, z] = meshgrid(linspace(-Lx, Lx, num_x), ...
                     linspace(-Ly, Ly, num_y), ...
                     linspace(0, Lz, num_z));
% 烟囱口 下标
center = round(Ly/dy)+1;
height = round(h/dz)+1;


%% 计算浓度函数
% 初始化
C = zeros(size(x));
figure(); 

% 模拟时间步进
for t = 0:dt:T
    C_new = C;
    % 更新浓度分布，除了边界
    for k = 2:num_z-1
        for j = 2:num_y-1
            for i = 2:num_x-1
                C_new(i,j,k) = C(i,j,k) + dt * D * ((C(i+1,j,k) - 2*C(i,j,k) + C(i-1,j,k)) / dx^2 ...
                    + (C(i,j+1,k) - 2*C(i,j,k) + C(i,j-1,k)) / dy^2 ...
                    + (C(i,j,k+1) - 2*C(i,j,k) + C(i,j,k-1)) / dz^2) ...
                    - dt * v * (C(i+1,j,k) - C(i-1,j,k)) / (2*dx) ...
                    - dt * w * (C(i,j,k+1) - C(i,j,k-1)) / (2*dz);
                if (C_new(i,j,k)<0)
                    C_new(i,j,k) = 0;
                end
            end
        end
    end
    
    %% 边界条件
    C_new(center, center, height) = C_new(center, center, height) + dt * C0;

    % 应用地面反射边界条件
    C_new(:, :, 1) = alpha * C(:, :, 2);

    %% 监测变量
    C_near = C_new(center-1:center+1, center-1:center+1, height);
    total_delta_kg = sum(C_new-C,"all");
    max_delta_kg = max(C_new-C, [], "all");
    
    % 更新浓度分布
    C = C_new;

    %% 绘图
    if (mod(t, 1)==0)
        % figure 1
        subplot(2,2,1)
        surf(y(:,:,1), x(:,:,1), C(:,:,1));% xy是反的，原因暂不明，fig2同
        colorbar;
        shading flat;

        xlabel('X/m');
        ylabel('Y/m');
        zlabel('C/(kg/m^3)');
        title(strcat('x-y plain at ground when t=', num2str(t)));

        xlim([-figy, figy]);
        ylim([-figx, figx]);

        % figure 2
        subplot(2,2,2)
        surf(y(:,:,height), x(:,:,height), C(:,:,height));
        colorbar;
        shading flat;

        xlabel('X/m');
        ylabel('Y/m');
        zlabel('C/(kg/m^3)');
        title(strcat('x-y plain at chimney when t=', num2str(t)));

        xlim([-figy, figy]);
        ylim([-figx, figx]);

        % figure 3
        subplot(2,2,3)
        surf(squeeze(y(:,center,:)), squeeze(z(:,center,:)), squeeze(C(:,center,:)));
        colorbar;
        shading flat;

        xlabel('X/m');
        ylabel('Z/m');
        zlabel('C/(kg/m^3)');
        title(strcat('x-z plain at chimney when t=', num2str(t)));

        xlim([-figy, figy]);
        ylim([0, figz]);

        % figure 4
        subplot(2,2,4)
        surf(squeeze(x(center,:,:)), squeeze(z(center,:,:)), squeeze(C(center,:,:)));
        colorbar;
        shading flat;

        xlabel('Y/m');
        ylabel('Z/m');
        zlabel('C/(kg/m^3)');
        title(strcat('y-z plain at chimney when t=', num2str(t)));

        xlim([-figx, figx]);
        ylim([0, figz]);

        drawnow;
    end
    
end

%% 完成模拟
fprintf('Simulation finished.\n');
