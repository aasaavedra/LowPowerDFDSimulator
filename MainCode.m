clear; clc; close all; % MATLAB Environment cleanup
tic; % Simulation runtime counter

%Direct-Fusion Drive Cold Plasma Dynamics and Propulsive Properties
%Simulator for Low-Power Experiments.

%% -----------------------------MAIN CODE---------------------------------
% ------------------------------------------------------------------------

%% I. INITIAL INFORMATION

% Physical data 
kB = 1.602e-19; % Boltzmann constant [eV/J]
Eion = 13.6; % Deuterium ionization energy [eV]
mD = 3.344e-27; % Mass of a Deuterium [kg]
mE = 9.10938356e-31; % Electron mass [kg]
L = 11.5; % System length (Gas box exit -> exhaust region) [m]
R0 = 0.05; % Radius at entrance [m]
mDot = 80e-6; % Injected mass flow rate[kg s^-1]
SHrat = 5/4; % Specific heat ratio of the plasma medium
Pw = 150e3; % Input power [W]
En = 0.1; % Neutral energy at the injection inlet [eV] 
g0 = 9.80665; %Acceleration to the Earth due to gravity [m/s^2]
alpha_inel=1; %Coefficient for excitation energy losses 

% Spatial domain
Nz = 3000; % Spatial grid points
dz = 1/(Nz-1); % Spatial time step
Z = linspace(0, 1, Nz); % Normalized spatial domain

% Error Parameters
errVar = 1; % Error quantifier
MaxErrorMain = 0.05; % Maximum error allowed
MaxError=MaxErrorMain;

% Normalization parameters
A_star = pi * R0^2; % Area normalization factor
T_star = Eion; % Temperature normalization factor
u_star = sqrt(kB * T_star / mD); % Ion velocity normalization factor
n_star = mDot / (A_star * mD * u_star); % Density normalization factor
t_star = L / u_star; % Time normalization factor
R_star = 1 / (t_star * n_star); % Ion. rate coefficient norm, factor
AQ_star = A_star * n_star * kB * T_star / t_star; % Plasma heat. norm. fac.
P_star = A_star * n_star * kB * T_star * u_star; %Power normalization fac.
F_star = mDot .* u_star; %Normalization of force magnitude;

% Computation of Area and Sigma Function
[Area, sig] = Asig(Nz, A_star, L);

% Energy source function, AQ.
zq = 7/23; lq = 3/23; % Plasma Heating param. (Gaussian profile)
AQnorm = exp(-((Z - zq) ./ lq).^2); % Norm. exp. for the plasma heat source
Pn = trapz(Z, AQnorm); % Integration of Q to obtain net input power
AQ0 = (Pw./P_star) / Pn; % Estimated amplitude of plasma heat source term
AQ = AQ0 .* AQnorm; % Resultant plasma heat source expression

% Ionization rate coefficient, Rion.
Rion0 = 1e-17 / Eion.^2 * sqrt(8 * kB * Eion / pi / mE);

% Precompute fRion values
Te_range = linspace(0.001, 10, 500); % Range of e temp. for precomputation
fRion_values = (Rion0./R_star)*arrayfun(@(Te) fRion(Te), Te_range);

% Simulation properties
CFL_lim = 0.95; % Courant–Friedrichs–Lewy Condition (CFL) limit
t = 0; % Initial simulation time
itr = 0; %Initial iteration count
pltItr = 1000; % Plot output every # of iterations
plt = pltItr; % Plot of initial snapshots

%% PLASMA EQUATIONS SOLVER: SETUP

%Initialization of dynamical variables
uzn = sqrt(2 * kB * En / mD) / u_star; %Initial neutral velocity (constant)
F = 0; % Thrust variable
eta = 0; % Propellant efficiency variable
etaU = 0; %Utilization Efficiency
Isp = 0; %Specific impulse
Ma = 0; %Mach Number
Elost = 0; %Energy lost to ionization

flag = 1; % Set '0' for default initial conditions or 
          %'1' for reuse of last simulation data
if flag == 1 && exist('LastInitData.mat', 'file')
    % Load previous simulation data
    load('LastInitData.mat', 'x1', 'x2', 'x3', 'x4', 'Zold');

    %Data interpolation to new grid, if necessary
    if length(Zold)~=length(Z)
    x1 = interp1(Zold, x1, Z, "spline", "extrap");
    x2 = interp1(Zold, x2, Z, "spline", "extrap");
    x3 = interp1(Zold, x3, Z, "spline", "extrap");
    x4 = interp1(Zold, x4, Z, "spline", "extrap");
    end

    %Redefine Max Error to avoid errors
    MaxErrorTemp = 0.000001;
    MaxError = MaxErrorTemp;

else %Use default initial conditions

    % Initial neutral density
    nn0_0 = 1 / uzn;
    nn0_1 = nn0_0 / 20;
    nn0 = nn0_0 * (1 - Z) + nn0_1 * Z;
        
    % Initial electron density
    ne0 = ones(1, Nz) * 0.1; 

    % Initial electron temperature
    Te0 = ones(1, Nz) * 0.1; 

    % Initial ion velocity parameters
    uzi0_0 = -0.25; % Initial constant
    uzi0_1 = 2; % Final constant
    uzi0 = uzi0_0 * (1 - Z) + uzi0_1 * Z; % Comp. of init. ion vel. values

    % Initialization of main physical variables
    x1 = Area .* ne0;
    x2 = Area .* ne0 .* uzi0; 
    x3 = 3/2 * Area .* ne0 .* Te0 + 1/2 * Area .* ne0 .* (uzi0 .^ 2);
    x4 = Area .* nn0;
end

% Computation of other physical variables
[ne, uzi, Te, nn, y1, y2, y3, Ape, ASp, cs] = NewVar(x1, x2, x3, x4, ...
    Area, SHrat, Te_range, fRion_values);

figure(2)
set(gcf, 'Position', get(0, 'Screensize'));

%% III. PLASMA EQUATIONS SOLVER: MAIN LOOP

% Main Loop
while errVar > MaxError
    % Plot of intermediate results
    if plt == pltItr
        figure(2)
        DataVis(Z, x1, x2, x3, x4, Area, SHrat, uzn, t, F, eta, Pw, ...
            mDot, En, Te_range, fRion_values, errVar, etaU, ...
            Isp, Ma, Elost);
        plt = 0;
    end

    %Restoration of original maximum error to the code
    if itr == 1000
        MaxError = MaxErrorMain;
    end

    % Time step update (based on CFL Criteron) 
    max_sp = max(abs(uzi) + cs); % Maximum possible speed
    dt = CFL_lim * dz / max_sp; % New time step

    % Flux computations (spatial derivatives)
    [k1x4, k1x1, k1x2, k1x3] = LLF(x4, x1, x2, x3, SHrat, dz, sig, ...
        Area, AQ, uzn, alpha_inel, Te_range, fRion_values);
    [k2x4, k2x1, k2x2, k2x3] = LLF(x4 + 0.5*dt*k1x4, x1 + 0.5*dt*k1x1, ...
        x2 + 0.5*dt*k1x2, x3 + 0.5*dt*k1x3, SHrat, dz, sig, Area, AQ, ...
        uzn,alpha_inel, Te_range, fRion_values);
    [k3x4, k3x1, k3x2, k3x3] = LLF(x4 + 0.5*dt*k2x4, x1 + 0.5*dt*k2x1, ...
        x2 + 0.5*dt*k2x2, x3 + 0.5*dt*k2x3, SHrat, dz, sig, Area, AQ, ...
        uzn, alpha_inel,Te_range, fRion_values);
    [k4x4, k4x1, k4x2, k4x3] = LLF(x4 + dt*k3x4, x1 + dt*k3x1, x2 ...
        + dt*k3x2, x3 + dt*k3x3, SHrat, dz, sig, Area, AQ, uzn, ...
        alpha_inel, Te_range, fRion_values);

    % Time Evolution
    % Temporal variables for main plasma variables' update
    dx1dt = (k1x1 + 2*k2x1 + 2*k3x1 + k4x1) / 6;
    dx2dt = (k1x2 + 2*k2x2 + 2*k3x2 + k4x2) / 6;
    dx3dt = (k1x3 + 2*k2x3 + 2*k3x3 + k4x3) / 6;
    dx4dt = (k1x4 + 2*k2x4 + 2*k3x4 + k4x4) / 6;

    x1 = x1 + dx1dt * dt;
    x2 = x2 + dx2dt * dt;
    x3 = x3 + dx3dt * dt;
    x4 = x4 + dx4dt * dt;
    
    % Boundary condition at left boundary (z=0)
    x1(1) = x1(2);
    x2(1) = -sqrt(2*SHrat/(SHrat+3)*x1(1)*x3(1));
    x3(1) = x3(2);
    x4(1) = (1 - x2(1)) / uzn;

    % Treatment of right boundary (z=0)
    x1(end) = x1(end-1);
    x2(end) = x2(end-1);
    x3(end) = x3(end-1);
    x4(end) = x4(end-1);
    
    % Temporal variables for dynamical properties' evolution
    [ne, uzi, Te, nn, y1, y2, y3, Ape, ASp, cs] = NewVar(x1, x2, x3,...
        x4, Area, SHrat, Te_range, fRion_values);

    % Time Update
    t = t + dt;
    plt = plt + 1;
    itr = itr + 1;
    clc; fprintf('Current time: %.4f s \n', t);
    fprintf('Total number of iterations: %.1f\n', itr);

    % Mass Flow data
    mFlow = y1 + x4 * uzn; % Mass flow variable
    errmFlow = mFlow / mFlow(1) - 1; % Mass flow relative error
    errVar = max(abs(errmFlow)); % Update of main error variable
    fprintf('Current error: %.4f %%\n', errVar*100);

    % Thrust computations
    F_i = y1(end) * uzi(end) * F_star; % Thrust from ions
    F_n = x4(end) * uzn^2 * F_star; % Thrust from neutrals
    F_e = 2/3 * (x3(end) - x2(end) * uzi(end) / 2) * F_star; % e-thrust
    F = F_i + F_e + F_n; % Net thrust

    %Specific Impulse
    Isp = F/(g0*mDot);
    
    % Propellant Efficiency
    eta = (F^2) / (2 * mDot * Pw);

    % Utilization Efficiency
    etaU = x2(end);

    %Mach Number
    Ma = uzi(end)/cs(end);
end

%Energy lost due to ionization
Elost = trapz(Z, ASp);
fprintf('\nEnergy spent for ionization: %.2f\n', Elost);

% Final plots
close all;
figure(3)
clf
set(gcf, 'Position', get(0, 'Screensize'));
DataVis(Z, x1, x2, x3, x4, Area, SHrat, uzn, t, F, eta, Pw, mDot, ...
    En, Te_range, fRion_values, errVar, etaU, Isp, Ma, Elost);

% Display of total simulation runtime
TotalTime = toc;
fprintf('Total simulation runtime: %.2f seconds\n', TotalTime);

% Save final simulation data for future use
Zold = Z; %Store Z value for interpolation to new grid
save('LastInitData.mat', 'x1', 'x2', 'x3', 'x4', 'Zold');

%Storage of resultant data
Zst=linspace(0, 1, 3000); %New standard grid for data analysis

%Interpolation of final values to new grid
Z=interp1(Z, Z, Zst, 'spline');
ne=interp1(Zold, ne, Zst, 'spline');
nn=interp1(Zold, nn, Zst, 'spline');
Te=interp1(Zold, Te, Zst, 'spline');
uzi=interp1(Zold, uzi, Zst, 'spline');
ASp=interp1(Zold, ASp, Zst, 'spline');
mFlow=interp1(Zold, mFlow, Zst, 'spline');

%Define matrix of final values.
OutputMatrix = [Z', ne', nn', Te', uzi', ASp', mFlow'];

%Save final data
writematrix(OutputMatrix, 'LatestResults.txt', 'Delimiter', 'tab');

%% -----------------------COMPLEMENTARY FUNCTIONS------------------------
% -----------------------------------------------------------------------

%% A. FUNCTION FOR AREA AND SIGMA PROCESSING.
function [Area, sig] = Asig(Nz, A_star, L)
    % Area Computation
    OldArea = load("DFD_Area.txt"); % Load area data of DFD system

    OldAreaLength = OldArea(:,1); % Store spatial grid of DFD system
    OldAreaLength = OldAreaLength / L; % Normalization of spatial grid

    OldAreaValues = OldArea(:,2); % Store area values of DFD system
    OldAreaValues = OldAreaValues / A_star; % Normalization of area values

    % Define new spatial grid for interpolation
    NewZ = linspace(min(OldAreaLength), max(OldAreaLength), Nz);
    % Interpolation of area data from old to new spatial grid
    Area = interp1(OldAreaLength, OldAreaValues, NewZ, 'spline');

    % Auxiliary Sigma Function Computation
    sigNum = [diff(log(Area))'; 0];
    sigDen = [diff(NewZ)'; 0];
    sig = sigNum ./ sigDen;

    % Output of area and sigma plots
    figure(1);
    set(gcf, 'Position', get(0, 'Screensize'));
    hold on;
    % Plot of area and sigma data, respectively
    plot(NewZ, Area, 'DisplayName', 'System Area');
    % Plot the Sigma Function data
    plot(NewZ, sig, 'DisplayName', 'Sigma Function');
    title('System Area and Sigma Function');
    xlabel('Normalized Length');
    ylabel('Area & Sigma Values');
    legend('show');
    hold off;
end

%% B. FUNCTION FOR COMPUTING DYNAMICAL QUANTITIES
function [ne, uzi, Te, nn, y1, y2, y3, Ape, ASp, cs] = NewVar(x1, x2, x3...
    , x4, Area, SHrat, Te_range, fRion_values)
    % Electron density (+)
    ne = max(x1, 0 .* x1) ./ Area;
    % Ion Velocity
    uzi = x2 ./ x1;
    % Area-scaled pressure term (+)
    Ape = max(2/3 * (x3 - x2 .* uzi / 2), 0 .* x1);
    % Electron temperature
    Te = Ape ./ x1; 
    % Neutral density
    nn = max(x4, 0 .* x4) ./ Area; 
    % Sonic speed
    cs = sqrt(abs(SHrat .* Te)); 
    % Ionization rate coefficient
    Rion = interp1(Te_range, fRion_values, Te, 'linear', 'extrap'); 
    % Area-scaled plasma source term
    ASp = max(x1.*x4.*Rion, 0.*x1) ./ Area;

    % Convected Variables
    y1 = x2; y2 = x2 .* uzi + Ape;
    y3 = SHrat / (SHrat - 1) * Ape .* uzi + x2 .* uzi .* uzi / 2;
end

%% C. FUNCTION FOR COMPUTING THE IONIZATION RATE COEFFICIENT
function [Rion] = fRion(T)
    N = 500; % Grid size
    % Fitting Parameters for ionization rate coefficient.
    A0 = 0.18450; A1 = -0.032226; A2 = -0.034539;
    A3 = 1.4003; A4 = -2.8115; A5 = 2.2986;
    Rion = 0 .* T; % Initialize ionization rate coefficient vector array

    for i = 1:length(T)
        T_i = T(i);
        % Integrand limits
        y = linspace (0.001, 1, N);
        
        % Definition of integrand
        x = 1 - T_i * log(y);
        u = 1 - 1 ./ x;
        inte = A0 * log(x) + u .* (A1 + u .* (A2 + u .* (A3 + u .* (A4...
            + u .* A5))));
    
        % Computation of integrand (trapezoidal method)
        Rion(i) = trapz(y, inte) / (sqrt(T_i) * exp(1 / T_i));
    end
end

%% D. DATA VISUALIZATION

function DataVis(Z, x1, x2, x3, x4, Area, SHrat, uzn, t, F, eta, Pw, ...
    mDot, En, Te_range, fRion_values, errVar, etaU, Isp, Ma, Elost)
    [ne, uzi, Te, nn, y1, ~, ~, ~, ASp, ~] = NewVar(x1, x2, x3, x4, ...
        Area, SHrat, Te_range, fRion_values);
    
    % Electron density
    subplot(2, 3, 1);
    hold on;
    plot(Z, ne, 'DisplayName', sprintf('t = %.4f', t));
    title('Electron Density');
    xlabel('Normalized Length');
    ylabel('Normalized Electron Density');
    %legend show;
    hold off;

    % Ion velocity
    subplot(2, 3, 2);
    hold on;
    plot(Z, uzi, 'DisplayName', sprintf('t = %.4f', t));
    title('Ion Velocity');
    xlabel('Normalized Length');
    ylabel('Normalized Ion Velocity');
    %legend show;
    hold off;

    % Electron temperature
    subplot(2, 3, 3);
    hold on;
    plot(Z, Te, 'DisplayName', sprintf('t = %.4f', t));
    title('Electron Temperature');
    xlabel('Normalized Length');
    ylabel('Normalized Electron Temperature');
    %legend show;
    hold off;

    % Neutral density
    subplot(2, 3, 4);
    hold on;
    plot(Z, nn, 'DisplayName', sprintf('t = %.4f', t));
    title('Neutral Density');
    xlabel('Normalized Length');
    ylabel('Neutral Density');
    %legend show;
    hold off;

    % Area-scaled plasma source
    subplot(2, 3, 5);
    hold on;
    plot(Z, ASp, 'DisplayName', sprintf('t = %.4f', t));
    title('Area-scaled Plasma Source');
    xlabel('Normalized Length');
    ylabel('Normalized Area-scaled P.S.');
    %legend show;
    hold off;
        
    % Mass Flow
    subplot(2, 3, 6);
    hold on;
    plot(Z, y1 + x4 .* uzn, 'DisplayName', sprintf('t = %.4f', t));
    title('Mass flow');
    xlabel('Normalized Length');
    ylabel('Normalized Mass Flow');
    %legend show;
    hold off;

    sgtitle(sprintf("- Plasma Dynamics of a Direct-Fusion Drive" + ...
        " -\nPw = %.2e W | mdot = %.0e kg/s | En = %.1f eV | SHrat =" + ...
        " %.2f | F = %.2f N | Elost = %.2f eV\n | eta = %.2f %%" + ...
        " | Isp = %.2f s | Ma = %.2f | etaU = %.2f %% | Error" + ...
        " = %.2f %% | t = %.2f s", Pw, mDot, En, SHrat, F, Elost, ...
        eta*100, Isp, Ma, etaU*100, errVar*100, t));

    
end

%% E. LOCAL LAX-FRIEDRICH SOLVER

function [kx4, kx1, kx2, kx3] = LLF(x4, x1, x2, x3, SHrat, dz, sig, ...
    Area, AQ, uzn, alpha_inel, Te_range, fRion_values)

    % A1. COMPUTATION OF TEMPORAL VARIABLES
    
    % Computation of spatial grid size
    N = length(x1);
    [kx1, kx2, kx3, kx4] = deal(zeros(1, N));

    % Computation of uzi, Ape and ASp values
    [~, uzi, ~, ~, y1, y2, y3, Ape, ASp, cs] = NewVar(x1, x2, x3, x4, ...
        Area, SHrat, Te_range, fRion_values);
    
    % LEFT-BOUNDARY TREATMENT : constant derivatives
    kx1(1) = 0; 
    kx2(1) = 0; 
    kx3(1) = 0; 
    kx4(1) = 0;
    
    % A2. NUMERICAL FLUX COMPUTATIONS WITH LLF SCHEME
    for i = 2:N-1
        % Computation of local wave speeds
        w_max_i = max(abs(uzi(i) + cs(i)), abs(uzi(i) - cs(i)));
        w_max_i1 = max(abs(uzi(i-1) + cs(i-1)), abs(uzi(i-1) - cs(i-1)));
        ai1 = max(w_max_i, w_max_i1);
    
        % Flux calculations with LLF scheme
        F1_minus = 0.5 * ((y1(i-1) + y1(i)) - ai1 * (x1(i) - x1(i-1)));
        F1_plus  = 0.5 * ((y1(i+1) + y1(i)) - ai1 * (x1(i+1) - x1(i)));
        F2_minus = 0.5 * ((y2(i-1) + y2(i)) - ai1 * (x2(i) - x2(i-1)));
        F2_plus  = 0.5 * ((y2(i+1) + y2(i)) - ai1 * (x2(i+1) - x2(i)));
        F3_minus = 0.5 * ((y3(i-1) + y3(i)) - ai1 * (x3(i) - x3(i-1)));
        F3_plus  = 0.5 * ((y3(i+1) + y3(i)) - ai1 * (x3(i+1) - x3(i)));
   
        % Derivative computation with diffusion terms
        kx1(i) = -(F1_plus - F1_minus) / dz + ASp(i);
        kx2(i) = -(F2_plus - F2_minus) / dz + Ape(i) * sig(i) ...
            + ASp(i) * uzn;
        kx3(i) = -(F3_plus - F3_minus) / dz + AQ(i) - ASp(i)...
            .*(alpha_inel-0.5.*uzn.^2);
        kx4(i) = (x4(i-1) - x4(i)) * uzn / dz - ASp(i);
    end
    
    % A3. RIGHT-BOUNDARY TREATMENT : constant derivatives
    kx1(N) = kx1(N-1);
    kx2(N) = kx2(N-1);
    kx3(N) = kx3(N-1);
    kx4(N) = kx4(N-1);
end


