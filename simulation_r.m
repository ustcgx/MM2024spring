% relations between the reflection coefficient of ground and the time spent for smoke to dissipate

clear;
clc;

D = 1;
v = 1;
h = 20;
T0 = 50;
C_thd = 50 * 1e-6 * 300 * 400 * 40;    % standard for average concentration of smoke: 50 mg/m^3
rr = [0.0, 0.2, 0.4, 0.6, 0.8];
T = zeros(size(rr));
C_sum = {}

%% simulation

for i = 1 : length(rr)
    r = rr(i);
    [T(i), C_sum{i}] = diffusion_reflection(r, T0, C_thd, D, v, h, false);
end

%% plot
% show the time spent for smoke to dissipate under different reflection coefficient of ground r

figure;
hold on;
for i = 1:length(C_sum)
    plot(C_sum{i} / 4.8);
end
hold off;

legend("r=0.0", "r=0.2", "r=0.4", "r=0.6", "r=0.8")
xlabel('Simulation time (s)')
ylabel('Average concentration in space (mg/(m^3))')

plot(rr, T);

xlabel('the reflection coefficient of ground')
ylabel('the time spent for smoke to dissipate (s), T0=50s')