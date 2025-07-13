set(0,'defaultaxesfontsize',20)
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultLineLineWidth',2);
set(0,'Defaultaxeslinewidth',2)

set(0,'DefaultFigureWindowStyle','docked')

c_c = 299792458;                % m/s TWM speed of light
c_eps_0 = 8.8542149e-12;        % F/m vaccum permittivity
c_eps_0_cm = c_eps_0/100;       % F/cm vaccum permittivity
c_mu_0 = 1/c_eps_0/c_c^2;       % Permiability of free space
c_q = 1.60217653e-19;           % Charge of an electon
c_hb = 1.05457266913e-34;       % Dirac / Reduced Planck constant
c_h = c_hb*2*pi;                % Planck constant

beta_r = 0;                     % De-tuning constant
beta_i = 0;                     % Gain Constant


beta_spe = .3e-5;
gamma = 1.0;
SPE = 7;

kappa0 = 0;                   % Coupling coefficient
kappaStart = 1/3;               % Constant defines starting position where coupling begins.
kappaStop = 2/3;                % Constant defines ending position where coupling stops.

InputParasL.E0 = 1e5;           % Amplitude of the input E-field / E_f
InputParasL.we = 0;             % Frequency of the complex sinusoidal modulation on the gaussian pulse
InputParasL.t0 = 200e-12;         % The constant we are shifting the time by
InputParasL.wg = 10e-13;         % Width of the Gaussian distribution
InputParasL.phi = 0;            % Initial Phase of the E_f / input E-field
InputParasR = 0;                % Placeholder variable for reverse propagation
InputParasL.rep = 500e-12;

n_g = 3.5;                      % Constant to control group velocity
vg = c_c/n_g *1e2;              % TWM cm/s group velocity

Lambda = 1550e-9;               % Wavelength of light
f0 = c_c/Lambda;

plotN = 10;                     % Divisior constant

% L = 1000e-6*1e2;                % length of the waveguide in cm
L = 1000e-6*1e2;                % length of the waveguide in cm


% Material Polarization Information
g_fwhm = 3.5e+012/10;           % Frequency
LGamma = g_fwhm*2*pi;
Lw0 = 0;
LGain = 0.015;                   % Gain Constant

XL = [0,L];                     % Start and End of the x-axis
%YL = [0,InputParasL.E0];       % Start and End of the y-axis
YL = [-InputParasL.E0,InputParasL.E0];       % Start and End of the y-axis

% Nz = 100;                       % Number of divisions
Nz = 100;                       % Number of divisions
dz = L/(Nz-1);                  % Distance between every point
dt = dz/vg;                     % Time taken to plot every point
fsync = dt*vg/dz;               % Equals 1, allows the Gaussian to be stable

Nt = floor(400*Nz);               % Time steps
tmax = Nt*dt;                   % Maximum time for simulation
t_L = dt*Nz;                    % time to travel length
z = linspace(0,L,Nz);           % Nz points, Nz-1 segments
time = nan(1,Nt);               % Time matrix with 1 row and Nt columns / row vector of Nt elements

InputL = nan(1,Nt);             % Matrix with 1 row and Nt columns / row vector of Nt elements
InputR = nan(1, Nt);            % Matrix with 1 row and Nt columns / row vector of Nt elements
OutputL = nan(1,Nt);            % Matrix with 1 row and Nt columns / row vector of Nt elements
OutputR = nan(1,Nt);            % Matrix with 1 row and Nt columns / row vector of Nt elements
Ef = zeros(size(z));            % Matrix with the same dimensions as z, all elements initialized to 0
Er = zeros(size(z));            % Matrix with the same dimensions as z, all elements initialized to 0

Ef1 = @SourceFct;               % Reference to SourceFct
ErN = @SourceFct;               % Reference to SourceFct

t = 0;                          % Set t to a starting value of 0
time(1) = t;                    % Sets the first element of the time vector to 0

InputL(1) = Ef1(t, InputParasL);  % Set initial value of InputL using the source function
InputR(1) = ErN(t, InputParasR);  % Set initial value of InputR using the source function

OutputR(1) = Ef(Nz);            % The end of the waveguide is the first value of the reflection (Right to Left)
OutputL(1) = Er(1);             % The end of the waveguide is the first value of the reflection (Left to Right)

Ef(1) = InputL(1);              % Initializes forward field at z = 0 (Input signal from the left)
Er(Nz) = InputR(1);             % Initializes backward field at z = L (Input signal from the right)

kappa = kappa0*ones(size(z));   % Creates an array of size z, where all indexes hold a value of kappa0
kappa(z<L*kappaStart) = 0;      % Sets the limit such that kappa is set to zero outside the interaction region
kappa(z>L*kappaStop) = 0;       % Sets the limit such that kappa is set to zero outside the interaction region.

Pf = zeros(size(z));            % Variable for the polarization of the material on the forward field
Pr = zeros(size(z));            % Variable for the polarization of the material on the reverse field

% Variables to hold field and polarization information
Efp = Ef;
Erp = Er;
Pfp = Pf;
Prp = Pr;

Nave = nan(1,Nt);
Ntr = 1e18;
N = ones(size(z))*Ntr;
Nave(1) = mean(N);

gain = vg*2.5e-16;
eVol = 1.5e-10*c_q;
Ion = 0.01e-9;
% Ion = 0.25e-9;
Ioff = 3e-9;
I_off = 0.024;
I_on = 0.1;
taun = 1e-9;
Zg = sqrt(c_mu_0/c_eps_0)/n_g;
EtoP = 1/(Zg*f0*vg*1e-2*c_hb);
alpha = 0;

figure('name', 'Fields')

% Forward field E_f
subplot(3,2,1)
plot(z*10000, real(Ef), 'r'); hold on
plot(z*10000, imag(Ef), 'r--');
plot(z*10000, real(Er), 'b--');
plot(z*10000, imag(Er), 'b');
hold off
xlim(XL*1e4)
ylim auto
xlabel('z (\mum)')
ylabel('E_f (V/\mum)')
legend('\Re(E_f)', '\Im(E_f)', '\Re(E_r)', '\Im(E_r)')


% Carrier Density N
subplot(3,2,2)
plot(z * 10000, N, 'r'); % Primary plot
xlim(XL * 1e4)
ylim([0, 5 * Ntr])
xlabel('z (\mum)')
ylabel('N')

% % Add a second y-axis for S
% yyaxis right
% plot(z * 10000, S, 'b'); % Plot S on the right axis
% ylabel('S') % Change the label to the appropriate units


% Average Carrier Density Over Time
subplot(3,2,3)
plot(time * 1e12, Nave, 'b');
xlim([0, Nt * dt * 1e12])
ylim([0, 5 * Ntr])
xlabel('time(ps)')
ylabel('Nave')

subplot(3,2,5)
plot(time * 1e12, real(OutputR), 'g'); hold on
plot(time * 1e12, real(OutputL), 'm--');
xlim([0, Nt * dt * 1e12])
ylim auto
xlabel('time(ps)')
ylabel('E (V/um)')
legend('Right Output', 'Left Output', 'Location', 'east')
hold off

for i = 2:Nt                    % 2 to 1000 in steps of 1

    t = dt*(i-1);
    time(i) = t;

    RL = 0.5;                               % The left side reflection coefficient
    RR = 0.5;                               % The right side reflection coefficient

    % Input
    InputL(i) = Ef1(t,0);
    % InputL(i) = Ef1(t, InputParasL);  % At time t, we input a signal characterized by InputParasL from the left
    InputR(i) = ErN(t, 0);            % At time t, we input no signal from the right (since InputParasR = 0)

    % Reflection
    Ef(1) = InputL(i) + RL*Er(1);     % Boundary condition at z = 0 (left side);
    Er(Nz) = InputR(i) + RR*Ef(Nz);   % Boundary condition at z = L (right side);

    S = (abs(Ef).^2 + abs(Er).^2).*EtoP*1e-6;

    if t < Ion || t > Ioff
        I_injv = I_off;
    else
        I_injv = I_on;
    end

    Stim = gain.*(N - Ntr).*S;
    N = (N + dt*(I_injv/ eVol - Stim))./(1+ dt/taun);
    Nave(i) = mean(N);

    gain_z = gain.*(N - Ntr)./vg;  % Compute gain coefficient
    beta_i = (gain_z - alpha)./2;  % Compute imaginary part of propagation constant
    beta = ones(size(z)).*(beta_r + 1i * beta_i); % Complex propagation constant
    exp_det = exp(-1i * dz * beta); % Phase shift due to propagation over dz

    Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz-1) + 1i*dz*kappa(2:Nz).*Er(2:Nz); % Forward Field Propagation
    Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz) + 1i*dz*kappa(2:Nz).*Ef(2:Nz);   % Reverse Field Propagation

    % Boundary Conditions
    Pf(1) = 0;     % zero polarization at the left boundary
    Pf(Nz) = 0;    % zero polarization at the right boundary
    Pr(1) = 0;     % zero polarization at the right boundary
    Pr(Nz) = 0;    % zero polarization at the left boundary
    Cw0 = -LGamma + 1i * Lw0;         % Defines the complex response function of the material.

    % Dispersion Calculations
    Tf = LGamma * Ef(1:Nz-2) + Cw0 * Pfp(2:Nz-1) + LGamma * Efp(1:Nz-2);  % Computes the forward polarization response based on previous field values.
    Pf(2:Nz-1) = (Pfp(2:Nz-1) + 0.5 * dt * Tf) ./ (1 - 0.5 * dt * Cw0);   % Updates the forward polarization field for every time step.
    Tr = LGamma * Er(3:Nz) + Cw0 * Prp(2:Nz-1) + LGamma * Erp(3:Nz);      % Computes the reverse polarization response based on previous field values.
    Pr(2:Nz-1) = (Prp(2:Nz-1) + 0.5 * dt * Tr) ./ (1 - 0.5 * dt * Cw0);   % Updates the reverse polarization field for every time step.

    Ef(2:Nz-1) = Ef(2:Nz-1) - LGain * (Ef(2:Nz-1) - Pf(2:Nz-1));  % Adjusts the forward electric field
    Er(2:Nz-1) = Er(2:Nz-1) - LGain * (Er(2:Nz-1) - Pr(2:Nz-1));  % Adjusts the reverse electric field

    % Output
    OutputR(i) = Ef(Nz) * (1 - RR); % Right output at z = L
    OutputL(i) = Er(1) * (1 - RL);  % Left output at z = 0

    A = sqrt(gamma*beta_spe*c_hb*f0*L*1e-2/taun)/(2*Nz);
    if SPE > 0
        eTf = ((randn(Nz,1)+1i*randn(Nz,1))*A).';
        eTr = ((randn(Nz,1)+1i*randn(Nz,1))*A).';
    else
        eTf = ((ones(Nz,1))*A).';
        eTr = ((ones(Nz,1))*A).';
    end

    EsF = eTf*abs(SPE).*sqrt(N.*1e6);
    Esr = eTr*abs(SPE).*sqrt(N.*1e6);

    Ef = Ef + EsF;
    Er = Er + Esr;

    % % FFT data from the outputs
    fftOutput1 = fftshift(fft(OutputR)); % Get FFT data for OutputR
    fftOutput2 = fftshift(fft(OutputL)); % Get FFT data for OutputL
    fftInput1 = fftshift(fft(InputL)); % Get FFT data for OutputL
    omega = fftshift(wspace(time));


    if mod(i,2000) == 0       % Only executed when i is multiple of plotN

        % Forward field E_f
        subplot(3,2,1)
        plot(z*10000, real(Ef), 'r'); hold on
        plot(z*10000, imag(Ef), 'r--');
        plot(z*10000, real(Er), 'b--');
        plot(z*10000, imag(Er), 'b');
        hold off
        xlim(XL*1e4)
        ylim auto
        xlabel('z (\mum)')
        ylabel('E_f (V/\mum)')
        legend('\Re(E_f)', '\Im(E_f)', '\Re(E_r)', '\Im(E_r)')


        % Carrier Density N
        subplot(3,2,2)
        plot(z * 10000, N, 'r'); % Primary plot
        xlim(XL * 1e4)
        ylim([0, 5 * Ntr])
        xlabel('z (\mum)')
        ylabel('N')

        % Add a second y-axis for S
        yyaxis right
        plot(z * 10000, S, 'b'); % Plot S on the right axis
        ylabel('S') % Change the label to the appropriate units

        % Average Carrier Density Over Time
        subplot(3,2,3)
        plot(time * 1e12, Nave, 'b');
        xlim([0, Nt * dt * 1e12])
        ylim([0, 2.5 * Ntr])
        xlabel('time(ps)')
        ylabel('Nave')

        % Input and Output Fields over time
        subplot(3,2,5)
        plot(time * 1e12, real(InputL), 'r'); hold on
        plot(time * 1e12, real(OutputR), 'g');
        plot(time * 1e12, real(InputR), 'b');
        plot(time * 1e12, real(OutputL), 'm--');
        xlim([0, Nt * dt * 1e12])
        ylim auto
        xlabel('time(ps)')
        ylabel('O')
        legend( 'Right Output', 'Left Output', 'Location', 'east')
        hold off


        % Output Field Spectrum (Magnitude in dB)
        subplot(3,2,4)
        plot(omega, 20*log10(abs(fftOutput1)), 'b'); hold on
        plot(omega, 20*log10(abs(fftInput1)), 'r'); %
        xlabel('Frequency (Hz)')
        ylabel('20 log_{10} |E|')
        legend('Output', 'Input')
        xlim([-0.05e14 0.05e14])
        hold off

        % Phase of Output Field
        subplot(3,2,6)
        phase1 = unwrap(angle(fftOutput1));  % Unwrap the phase1
        phase2 = unwrap(angle(fftInput1));  % Unwrap the phase for Input
        plot(omega, phase1);
        hold on;
        plot(omega, phase2);
        hold off;
        xlabel('Frequency (Hz)')
        ylabel('phase (E)')
        legend('fftOutput1', 'fftInput1');

        pause(0.01)


    end


    % Update Previous Values
    Efp = Ef;
    Erp = Er;
    Pfp = Pf;
    Prp = Pr;
end
% disp('Simulation finished. Running FFT analysis...');

% %% --- Compute FFT from outputR during transient and steady-state
% transient_range = 1:round(0.005*Nt);     
% steady_range = round(0.0001*Nt):round(1*Nt);      
% 
% output_transient = OutputR(transient_range);
% output_steady = OutputR(steady_range);
% 
% time_transient = time(transient_range);
% time_steady = time(steady_range);
% 
% omega_transient = fftshift(wspace(time_transient));
% omega_steady = fftshift(wspace(time_steady));
% 
% fft_transient = fftshift(fft(output_transient));
% fft_steady = fftshift(fft(output_steady));
% 
% % --- Side-by-Side Transient and Steady-State FFTs
% figure('Name', 'Transient vs Steady-State', 'Color', 'w');
% 
% % Transient FFT
% subplot(1,2,1);
% plot(omega_transient, 20*log10(abs(fft_transient)), 'b');
% xlabel('Frequency (Hz)');
% ylabel('20 log_{10} |E|');
% title('Modes during Transient ');
% grid on;
% xlim([-0.1e14 0.1e14]);
% 
% % Steady-State FFT
% subplot(1,2,2);
% plot(omega_steady, 20*log10(abs(fft_steady)), 'r');
% xlabel('Frequency (Hz)');
% ylabel('20 log_{10} |E|');
% title('Modes Steady-State');
% grid on;
% xlim([-0.02e14 0.02e14]);




