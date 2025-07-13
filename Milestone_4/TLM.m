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

kappa0 = 0;                   % Coupling coefficient
kappaStart = 1/3;               % Constant defines starting position where coupling begins.
kappaStop = 2/3;                % Constant defines ending position where coupling stops.

InputParasL.E0 = 1e5;           % Amplitude of the input E-field / E_f
InputParasL.we = 0;             % Frequency of the complex sinusoidal modulation on the gaussian pulse
InputParasL.t0 = 2e-12;         % The constant we are shifting the time by
InputParasL.wg = 5e-13;         % Width of the Gaussian distribution
InputParasL.phi = 0;            % Initial Phase of the E_f / input E-field
InputParasR = 0;                % Placeholder variable for reverse propagation

n_g = 3.5;                      % Constant to control group velocity
vg = c_c/n_g *1e2;              % TWM cm/s group velocity

Lambda = 1550e-9;               % Wavelength of light

plotN = 10;                     % Divisior constant

L = 1000e-6*1e2;                % length of the waveguide in cm

% Material Polarization Information
g_fwhm = 3.5e+012/10;           % Frequency
LGamma = g_fwhm*2*pi;            
Lw0 = 0;                        
LGain = 0.01;                   % Gain Constant

XL = [0,L];                     % Start and End of the x-axis
YL = [-InputParasL.E0,InputParasL.E0];       % Start and End of the y-axis

Nz = 500;                       % Number of divisions
dz = L/(Nz-1);                  % Distance between every point
dt = dz/vg;                     % Time taken to plot every point
fsync = dt*vg/dz;               % Equals 1, allows the Gaussian to be stable

Nt = floor(2*Nz);               % Time steps
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

figure('name', 'Fields')

subplot(3,1,1)
plot(z*1000, real(Ef), 'r');
hold off
xlabel('z(\mum)')
ylabel('E_f')

subplot(3,1,2)
plot(z*1000, real(Er), 'b');
xlabel('z(\mum)')
ylabel('E_r')
hold off

subplot(3,1,3)
plot(time*1e12, real(InputL), 'r'); hold on
plot(time*1e12, real(OutputR), 'r--');
plot(time*1e12, real(InputR), 'b'); hold on
plot(time*1e12, real(OutputL), 'b--');
xlabel('time(ps)')
ylabel('E')
hold off

for i = 2:Nt                    % 2 to 1000 in steps of 1

    t = dt*(i-1);
    time(i) = t;

    RL = 0;                               % The left side reflection coefficient
    RR = 0;                               % The right side reflection coefficient

    beta = ones(size(z))*(beta_r+1i*beta_i); % Complex propagation constant
    exp_det = exp(-1i*dz*beta);              % Phase shift due to propagation over a distance dz

    % Input
    InputL(i) = Ef1(t, InputParasL);  % At time t, we input a signal characterized by InputParasL from the left
    InputR(i) = ErN(t, 0);            % At time t, we input no signal from the right (since InputParasR = 0)

    % Reflection
    Ef(1) = InputL(i) + RL*Er(1);     % Boundary condition at z = 0 (left side);
    Er(Nz) = InputR(i) + RR*Ef(Nz);   % Boundary condition at z = L (right side);

    Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz-1) + 1i*dz*kappa(2:Nz).*Er(2:Nz); % Forward Field Propagation
    Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz) + 1i*dz*kappa(2:Nz).*Ef(2:Nz);   % Reverse Field Propagation

    % Boundary Conditions
    Pf(1) = 0;     % zero polarization at the left boundary                    
    Pf(Nz) = 0;    % zero polarization at the right boundary                
    Pr(1) = 0;     % zero polarization at the right boundary                  
    Pr(Nz) = 0;    % zero polarization at the left boundary                    
    Cw0 = -LGamma + 1i * Lw0;         % Defines the complex response function of the material.

    % Dispersion Calculations
    % Backward Euler Polarization Update
    % Forward polarization
    Tf = LGamma * Efp(2:Nz-1);  % Source term for forward polarization
    Pf(2:Nz-1) = (Pfp(2:Nz-1) + dt * Tf) ./ (1 - dt * Cw0);  % Forward polarization update

    % Backward polarization
    Tr = LGamma * Erp(2:Nz-1);  % Source term for backward polarization
    Pr(2:Nz-1) = (Prp(2:Nz-1) + dt * Tr) ./ (1 - dt * Cw0);  % Backward polarization update 
 
    Ef(2:Nz-1) = Ef(2:Nz-1) - LGain * (Ef(2:Nz-1) - Pf(2:Nz-1));  % Adjusts the forward electric field 
    Er(2:Nz-1) = Er(2:Nz-1) - LGain * (Er(2:Nz-1) - Pr(2:Nz-1));  % Adjusts the reverse electric field 
    
    % Output
    OutputR(i) = Ef(Nz) * (1 - RR); % Right output at z = L
    OutputL(i) = Er(1) * (1 - RL);  % Left output at z = 0

    % FFT data from the outputs
    fftOutput1 = fftshift(fft(OutputR)); % Get FFT data for OutputR
    fftOutput2 = fftshift(fft(OutputL)); % Get FFT data for OutputL
    fftInput1 = fftshift(fft(InputL)); % Get FFT data for OutputL
    omega = fftshift(wspace(time));

    if mod(i,plotN) == 0       % Only executed when i is multiple of plotN

        % Forward Propagation of the Gaussian Pulse
        subplot(3,2,1)
        plot(z*10000,real(Ef),'r'); hold on
        plot(z*10000,imag(Ef),'r--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_f (V/um)')
        legend('\Re','\Im')
        hold off

        % Reverse (Reflection) of the Gaussian Pulse
        subplot(3,2,3)
        plot(z*10000, real(Er), 'b'); hold on
        plot(z*10000, imag(Er), 'b--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_r (V/um)')
        legend('\Re', '\Im')
        hold off

        % Plot showing when the time when the input and output pulse were detected
        subplot(3,2,5);
        plot(time*1e12, real(InputL), 'r'); hold on
        plot(time*1e12, real(OutputR), 'g');
        plot(time*1e12, real(InputR), 'b');
        plot(time*1e12, real(OutputL), 'm');
        xlim([0, Nt*dt*1e12])
        ylim(YL)
        xlabel('time(ps)')
        ylabel('E (V/um)')
        legend('Left Input', 'Right Output', 'Right Input', 'Left Output' ...
            , 'Location', 'east')
        hold off

        % Plot showing the spectral content of OutputR and OutputL
        subplot(3,2,2);
        plot(omega, abs(fftOutput1));
        hold on; %
        plot(omega, abs(fftOutput2));
        plot(omega, abs(fftInput1));
        hold off; %
        xlabel('Frequency (THz)');
        ylabel('|E| (V/um)');
        xlim([-0.1e14, 0.1e14]);
        legend('fftOutputR','fftOutputL','fftInputL');

        % Plot showing the phase of OutputR and OutputL
        subplot(3,2,4);
        phase1 = unwrap(angle(fftOutput1));  % Unwrap the phase1
        phase2 = unwrap(angle(fftOutput2));  % Unwrap the phase2
        phase3 = unwrap(angle(fftInput1));  % Unwrap the phase for Input

        plot(omega, phase1);
        hold on;
        plot(omega, phase2);
        plot(omega, phase3);
        hold off;
        xlabel('Frequency (THz)');
        ylabel('Phase (E) (Rads)');
        xlim([-1.5e14, 1.5e14]);
        legend('fftOutputR','fftOutputL','fftInputL');

        pause(0.01)
    end

    % Update Previous Values
    Efp = Ef;
    Erp = Er;
    Pfp = Pf;
    Prp = Pr;
end