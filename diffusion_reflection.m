% MATLAB simulation program
% modified from diffusion.m to take ground reflection into consideration
% There are 5 main diffenrences:
% 1. add argument r, T0, C_thd, remove return value res
% 2. change argument T to a return value, it now means the time spent for smoke to dissipate
% 3. line 58 change for to while, simulate until the total concentration drops to the threshold
% 4. line 77 add if, limit the working time of chimeny
% 5. line 82 consider the influence of ground absorption in boundary conditions
function T = diffusion_reflection(r, T0, C_thd, D, v, h, show_fig, C_init, t_init)
    %% Parameter Description
    % r: the reflection coefficient of ground, 0<=r<=1.0
    % T0: the working time of the chimney
    % C_thd: the threshold of total smoke concentration in space

    %% Other parameter settings
    Lx = 300; % x-axis simulation range in meters
    Ly = 300; % y-axis simulation range in meters
    Lz = 100; % z-axis simulation range in meters

    figx = 200; % x-axis visualization range in meters
    figy = 200; % y-axis visualization range in meters
    figz = 100; % z-axis visualization range in meters

    dx = 1; % x-direction grid size in meters
    dy = 1; % y-direction grid size in meters
    dz = 1; % z-direction grid size in meters
    dt = 0.1; % Time step in seconds

    C0 = 1e-2; % Constant concentration emitted by the chimney in kg/s

    %% Create grid
    num_x = round(2*Lx/dx)+1;
    num_y = round(2*Ly/dy)+1;
    num_z = round(Lz/dz)+1;
    [x, y, z] = meshgrid(linspace(-Lx, Lx, num_x), ...
                        linspace(-Ly, Ly, num_y), ...
                        linspace(0, Lz, num_z));
    % Chimney outlet Subscript
    center = round(Ly/dy)+1;
    height = round(h/dz)+1;

    %% Calculate concentration function
    % Initialization
    if (nargin<=7)
        % C, t not initiated
        C = zeros(size(x));
        t = 0;
    else
        % initial C, t given
        C = C_init;
        t = t_init;
    end

    if (show_fig)
        figure(); 
    end
    % Simulation time stepping
    while (t <= T0 || sum(C, "all") > C_thd)
        t = t + dt;
        C_new = C;
        % Update concentration distribution, except for boundaries
        for k = 2:num_z-1
            for j = 2:num_y-1
                for i = 2:num_x-1
                    C_new(i,j,k) = C(i,j,k) + dt * D * ((C(i+1,j,k) - 2*C(i,j,k) + C(i-1,j,k)) / dx^2 ...
                        + (C(i,j+1,k) - 2*C(i,j,k) + C(i,j-1,k)) / dy^2 ...
                        + (C(i,j,k+1) - 2*C(i,j,k) + C(i,j,k-1)) / dz^2) ...
                        - dt * v * (C(i+1,j,k) - C(i-1,j,k)) / (2*dx);
                    if (C_new(i,j,k)<0)
                        C_new(i,j,k) = 0;
                    end
                end
            end
        end
    
        %% Boundary conditions
        if (t <= T0)
            C_new(center, center, height) = C_new(center, center, height) + dt * C0;
        end
        
        % Apply ground reflection boundary conditions (zero flux) assuming the concentration remains unchanged after reflection
        C_new(:, :, 1) = (1-r) * C(:, :, 2);    % consider the influence of ground absorption in boundary conditions

        %% Monitor variables
        % C_near = C_new(center-1:center+1, center-1:center+1, height)
        total_delta_kg = sum(C_new-C,"all");
        % max_delta_kg = max(C_new-C, [], "all")
        
        % Update concentration distribution
        C = C_new;

        %% Plotting
        if (show_fig && mod(t, 1)==0)
            % figure 1
            subplot(2,2,1);
            visualization(y(:,:,1), x(:,:,1), C(:,:,1), 'X/m', 'Y/m', 'C/(kg/m^3)', ...
                [-figy figy], [-figx figx], ...
                strcat('x-y plain at ground when t=', num2str(t)), ...
                [0, 3e-5]);
            
            % figure 2
            subplot(2,2,2);
            visualization(y(:,:,height), x(:,:,height), C(:,:,height), 'X/m', 'Y/m', 'C/(kg/m^3)', ...
                [-figy figy], [-figx figx], ...
                strcat('x-y plain at chimney when t=', num2str(t)));

            % figure 3
            subplot(2,2,3);
            visualization(squeeze(y(:,center,:)), squeeze(z(:,center,:)), squeeze(C(:,center,:)), 'X/m', 'Z/m', 'C/(kg/m^3)', ...
                [-figy figy], [0, figz], ...
                strcat('x-z plain at chimney when t=', num2str(t)));

            % figure 4
            subplot(2,2,4);
            visualization(squeeze(x(center,:,:)), squeeze(z(center,:,:)), squeeze(C(center,:,:)), 'Y/m', 'Z/m', 'C/(kg/m^3)', ...
                [-figx figx], [0, figz], ...
                strcat('y-z plain at chimney when t=', num2str(t)));

            %% Save
            % if (mod(t, 50)==0 && t>0)
            %     savefig(strcat('t=', num2str(t),'.fig'))
            % end
        end

        %% Convergence
        if (total_delta_kg<1e-6)
            break;
        end
    end

    %% Simulation finished
    % if (show_fig)
    %     % for specific task, please modify the title
    %     savefig(strcat('t=', num2str(t),'.fig'))
    % end
    fprintf('Simulation finished.\n');
    T = t - T0;
end