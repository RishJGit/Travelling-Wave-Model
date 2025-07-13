set(0,'defaultaxesfontsize',20)
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultLineLineWidth',2);
set(0,'Defaultaxeslinewidth',2)

set(0,'DefaultFigureWindowStyle','docked')

c_c = 299792458;                % m/s TWM speed of light
c_eps_0 = 8.8542149e-12;        % F/m vaccum permittivity
c_eps_0_cm = c_eps_0/100;       % F/cm vaccum permittivity
c_mu_0 = 1/c_eps_0/c_c^2;       % Permeability of free space
c_q = 1.60217653e-19;           % Charge of an electon
c_hb = 1.05457266913e-34;       % Dirac / Reduced Planck constant
c_h = c_hb*2*pi;                % Planck constant

beta_r = 80;                    % Pulse De-tuning Constant
beta_i = 8;                     % Pulse Gain Constant

InputParasL.E0 = 1e5;           % Amplitude of the input E-field / E_f
InputParasL.we = 0;        % Frequency of the complex sinusoidal modulation on the gaussian pulse 
InputParasL.t0 = 2e-12;         % Time delay of the input pulse
InputParasL.wg = 5e-13;         % Width of the Gaussian distribution
InputParasL.phi = 0;           % Initial Phase of the E_f / input E-field
InputParasR = 0;                % Placeholder variable for reverse propagation

n_g = 3.5;                     % Constant to control group velocity
vg = c_c/n_g *1e2;             % TWM cm/s group velocity
Lambda = 1550e-9;              % Wavelength of light

plotN = 10;                    % Divisor constant

L = 1000e-6*1e2;               % length of the waveguide in cm
XL = [0,L];                    % Start and End of the x-axis
YL = [-2*InputParasL.E0,2*InputParasL.E0];       % Start and End of the y-axis

Nz = 500;                      % Number of divisions 
dz = L/(Nz-1);                 % Spatial step size in the waveguide (cm)
dt = dz/vg;                    % Time step size for wave propagation (s)
fsync = dt*vg/dz;              % Stability factor, ensures synchronization

Nt = floor(2*Nz);              % Maximum Number of time steps to simulate up to
tmax = Nt*dt;                  % Maximum time (in s) to run the simulation up to
t_L = dt*Nz;                   % time to travel length

z = linspace(0,L,Nz).';        % z goes from 0 to Nz points, in L segments
time = nan(1,Nt);              % Time matrix with 1 row and Nt columns / row vector of Nt elements
InputL = nan(1,Nt);            % Matrix with 1 row and Nt columns / row vector of Nt elements
InputR = nan(1, Nt);           % Matrix with 1 row and Nt columns / row vector of Nt elements
OutputL = nan(1,Nt);           % Matrix with 1 row and Nt columns / row vector of Nt elements
OutputR = nan(1,Nt);           % Matrix with 1 row and Nt columns / row vector of Nt elements

Ef = zeros(size(z));          % Forward Electric field at each of spatial location. Each Element initialized to 0.
Er = zeros(size(z));          % Reverse Electric field in at each of spatial location. Each Element initialized to 0.

Ef1 = @SourceFct;              % Function handle for source function
ErN = @SourceFct;              % Function handle for source function

t = 0;                         % Set t to a starting value of 0
time(1) = t;                   % Sets the first element of the time vector to 0

InputL(1) = Ef1(t, InputParasL);  % Set initial value of InputL using the source function
InputR(1) = ErN(t, InputParasR);  % Set initial value of InputR using the source function

OutputR(1) = Ef(Nz);           % The end of the waveguide is the first value of the reflection (Right to Left)
OutputL(1) = Er(1);            % The end of the waveguide is the first value of the reflection (Left to Right)

Ef(1) = InputL(1);
Er(Nz) = InputR(1);

figure('name', 'Fields')
subplot(3,1,1)
plot(z*1000, real(Ef), 'r');
hold off
xlabel('z(\mum)')
ylabel('E_f (V/um)')
subplot(3,1,2)
plot(z*1000, real(Er), 'b');
xlabel('z(\mum)')
ylabel('E_r (V/um)')
hold off
subplot(3,1,3)
plot(time*1e12, real(InputL), 'r'); hold on
plot(time*1e12, real(OutputR), 'r--');
plot(time*1e12, real(InputR), 'b'); hold on
plot(time*1e12, real(OutputL), 'b--');
xlabel('time(ps)')
ylabel('E')

hold off

for i = 2:Nt                % 2 to 1000 in steps of 1
    
    t = dt*(i-1);
    time(i) = t;

    RL = 0; % The reflection coefficient
    RR = 0.9i; % The reflection coefficient
    
    beta = ones(size(z))*(beta_r+1i*beta_i); %
    exp_det = exp(-1i*dz*beta);
    
    % Input
    InputL(i) = Ef1(t, InputParasL);
    InputR(i) = ErN(t, 0);
    
    % Reflection 
    Ef(1) = InputL(i) + RL*Er(1); % Left 
    Er(Nz) = InputR(i) + RR*Ef(Nz); % Right
    % Forward Propagation
    Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz-1);
    Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz);
    
    % Output
    OutputR(i) = Ef(Nz)*(1-RR);  
    OutputL(i) = Er(1)*(1-RL);  

    fftOutput = fftshift(fft(OutputR));    % Shifting the output and saving it in as a vector)
    omega = fftshift(wspace(time));        % 

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
        
        % Plot the FFT (Gain vs Frequency) of the Gaussion (or Gaussian Modulated) Pulse
        subplot(3,2,2);
        plot(omega, abs(fftOutput));
        xlabel('Frequency (THz)');
        ylabel('|E| (V/um)');
        xlim([-1.5e14, 1.5e14]);

        % Plot the Phase of the Gaussion (or Gaussian Modulated) Pulse
        subplot(3,2,4);
        phase = unwrap(angle(fftOutput));  % Unwrap the phase
        plot(omega, phase);
        xlabel('Frequency (THz)');
        ylabel('Phase (E) (Rads)');
        xlim([-1.5e14, 1.5e14]);
        
        
        pause(0.01)
    end
end