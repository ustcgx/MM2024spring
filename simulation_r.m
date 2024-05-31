% relations between the reflection coefficient of ground and the time spent for smoke to dissipate

clear;
clc;

D = 1;
v = 1;
h = 20;
T0 = 50;    % the working time of the chimeny
C_thd = 50 * 1e-9 * 300 * 600 * 40;    % standard for average concentration of smoke: 50 ug/m^3
rr = [0.0, 0.2, 0.4, 0.6, 0.8];
T = zeros(size(rr));
C_sum = {};

%% simulation

for i = 1 : length(rr)
    r = rr(i);
    [T(i), C_sum{i}] = diffusion_reflection(r, T0, C_thd, D, v, h, true);
end

%% plot
% show the time spent for smoke to dissipate under different reflection coefficient of ground r

figure;
hold on;
for i = 1:length(C_sum)
    % change kg to ug/m^3
    plot(C_sum{i} * 1000 / 7.2, 'DisplayName', ['r=', num2str(rr(i))]);
end
hold off;

legend();
xlabel('Simulation time (s)')
ylabel('Average concentration in space (ug/(m^3))')

figure;
T = round(T);
plot(rr, T, '-o');
text(rr, T, num2str(T'), 'VerticalAlignment','bottom','HorizontalAlignment','right');

xlabel('the reflection coefficient of ground')
ylabel('the time spent for smoke to dissipate (s), T0=50s')