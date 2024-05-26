% relations between chimney height and max concentration at ground

clear;
clc;

D = 1;
v = 1;
TT = [50, 100, 200, 400];
hh = [10, 20, 30, 40, 50];
max_c = zeros(4, 5);

%% simulation

for i = 1 : length(hh)
    h = hh(i);
    for j = 1 : length(TT)
        T = TT(j);
        if (j==1)
            C = diffusion(D, v, h, T, false);
        else
            C = diffusion(D, v, h, T, false, C_last, T_last);
        end
        T_last = T;
        C_last = C;
        max_c(j, i) = max(C_last(:,:,1), [], "all");
    end
end

%% plot
% show max concentration at ground at diffenrent t under different chimney height t

figure;
plot(hh, max_c);

legend('t=50', 't=100', 't=200', 't=400')
xlabel('Chimney height (m)')
ylabel('Max concentration at ground (kg/(m^3·s))')
