set(0, 'defaultaxesfontsize', 20)
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesLineWidth', 2)
set(0, 'DefaultFigureWindowStyle', 'docked')

c_c = 299792458;
c_eps_0 = 8.8542149e-12;
c_eps_0_cm = c_eps_0/100;
c_mu_0 = 1 / c_eps_0 / c_c^2;
c_hb = 1.05457266913e-34;
c_h = c_hb * 2 * pi;

InputParasL.E0 = 1e5;
InputParasL.we = 0;
InputParasL.t0 = 2e-12;
InputParasL.wg = 2e-13;
InputParasL.phi = 0;
InputParasR = 0;

RL = 0;
RR = 0;

n_g = 3.5;
vg = c_c / n_g * 1e2;
Lambda = 1550e-9;

plotN =10;

L = 1000e-6 * 1e2;
XL = [0,L];
YL = [-2.5InputParasL.E0,2.5InputParasL.E0];

Nz = 500;
dz = L / (Nz - 1);
dt = dz / vg;
fsync = dt * vg / dz;

Nt = floor(2 * Nz);
tmax = Ntdt;
t_L = dtNz;

z = linspace(0, L, Nz).';

kappa0 = 100;
central_freq = L / 2;
bandwidth = L / 10;
kappa = kappa0 * exp(-((z - central_freq) / bandwidth).^2) .* cos(2 * pi * (z / L) * 10);
kappa(z < L/3) = 0;
kappa(z > 2L/3) = 0;

figure('Name', 'Grating Coupling Coefficient');
plot(z 1e4, kappa, 'b', 'LineWidth', 2);
xlabel('z (μm)');
ylabel('κ');
title('Grating Coupling Coefficient κ as a Function of z');
grid on;

time = nan(1, Nt);
InputL = nan(1, Nt);
InputR = nan(1, Nt);
OutputL = nan(1, Nt);
OutputR = nan(1, Nt);
Ef = zeros(size(z));
Er = zeros(size(z));

Ef1 = @SourceFct;
ErN = @SourceFct;
t = 0;
time(1) = t;
InputL(1) = Ef1(t, InputParasL);
InputR(1) = ErN(t, InputParasR);
OutputR(1) = Ef(Nz);
OutputL(1) = Er(1);
Ef(1) = InputL(1);
Er(Nz) = InputR(1);

beta_r = 0;
beta_i = 0;

beta = ones(size(z)) * (beta_r + 1i * beta_i);
exp_det = exp(-1i * dz * beta);

for i = 2:Nt
    t = dt * (i - 1);
    time(i) = t;
    InputL(i) = Ef1(t, InputParasL);
    InputR(i) = ErN(t, InputParasR);
    Ef(1) = InputL(i) + RL * Er(1);
    Er(Nz) = InputR(i) + RR * Ef(Nz);
    Ef_prev = Ef;
    Ef(2:Nz) = fsync * exp_det(1:Nz-1).Ef(1:Nz-1) + 1i dz * kappa(2:Nz) .* Er(2:Nz);
    Er(1:Nz-1) = fsync * exp_det(2:Nz).Er(2:Nz) + 1i dz * kappa(1:Nz-1) .* Ef_prev(1:Nz-1);
    OutputR(i) = Ef(Nz) * (1 - RR);
    OutputL(i) = Er(1) * (1 - RL);
end

fftInputL = fftshift(fft(InputL));
fftOutputR = fftshift(fft(OutputR));
fftOutputL = fftshift(fft(OutputL));
omega = fftshift(wspace(time));
omega_THz = omega / (2 * pi * 1e12);

figure('Name', 'FFT Magnitude Analysis');
plot(omega_THz, abs(fftInputL), 'r', 'DisplayName', 'Input');
hold on;
plot(omega_THz, abs(fftOutputR), 'g', 'DisplayName', 'Right Output');
plot(omega_THz, abs(fftOutputL), 'b', 'DisplayName', 'Left Output');
hold off;
xlabel('Frequency (THz)');
xlim([-5, 5]);
ylabel('Magnitude (V/m)');
title('FFT Magnitude of Input and Outputs');
legend;
grid on;

figure('Name', 'FFT Phase Analysis');
plot(omega_THz, unwrap(angle(fftInputL)), 'r', 'DisplayName', 'Input');
hold on;
plot(omega_THz, unwrap(angle(fftOutputR)), 'g', 'DisplayName', 'Right Output');
plot(omega_THz, unwrap(angle(fftOutputL)), 'b', 'DisplayName', 'Left Output');
hold off;
xlabel('Frequency (THz)');
ylabel('Phase (rad)');
title('FFT Phase of Input and Outputs');
legend;
grid on;
