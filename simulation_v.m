% relations between wind speed and max concentration at ground

clear;
clc;

D = 1;
h = 20;
TT = [30, 60, 90];
vv = [0.1, 0.5, 1, 2, 4];

%% simulation

for i = 1 : length(vv)
    v = vv(i);
    for j = 1 : length(TT)
        T = TT(j);
        if (j==1)
            C = diffusion(D, v, h, T, true);
        else
            C = diffusion(D, v, h, T, true, C_last, T_last);
        end
        T_last = T;
        C_last = C;
    end
end

% directly save figure at different t under different v.