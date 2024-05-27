% consider the influence of gravity, 
% call diffusion_gravity instead of diffusion !!!
 

clear;
clc;

D = 1;
h = 20;
v = 1;
TT = [30, 60, 90];
gg = [-1, 1]; % g<0 means that there is a wind towards positive z-axis

%% simulation

for i = 1 : length(gg)
    g = gg(i);
    for j = 1 : length(TT)
        T = TT(j);
        if (j==1)
            C = diffusion_gravity(g, D, v, h, T, true);
        else
            C = diffusion_gravity(g, D, v, h, T, true, C_last, T_last);
        end
        T_last = T;
        C_last = C;
    end
end

% when call diffusion_gravity, please save figure
% directly save figure at different t under different g.
% to show that g contributes to the concentration on the ground