%% Computation of the initial data
%% INPUTS:
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%%          Physics2D: Structure containing variables concerning the physics of the problem in 2D (structure) (see Physics2D_Var2d.m)
%% INPUTS(OPTIONAL):
%%          InitialData_Choice: Variable containing the choice between type of computations for the initial data (double, Default: 1)
%%          (Must either be: 1 to compute directly a centered gaussian, 2 to compute directly the Thomas-Fermi approximation)
%%          Rabi: radio-frequency field frequency
%%          delta: detuning value to select the initial spin polarization
%%          X0,Y0: Coordinates of the center of Gaussians or Thomas-Fermi approximation (vector or double, Default: 0)
%%          gamma_x,gamma_y:Parameters for the centered gaussian (double, Default: 1) (see GaussianInitialData2d.m)
%% OUTPUT:
%%          Phi_0: Initial data computed (cell array)
%% FUNCTIONS USED:
%%          GaussianInitialData2d: To compute the centered gaussian (line 59,76 and 107)
%%          Thomas_Fermi2d: To compute the Thomas-Fermi approximation (line 67 and 95)
%%          CNSP_CNGF2d: To compute initial data using the CNSP-CNFG method (line 101 and 113)

function [Phi_0] = InitialData_Var2d_custom(varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addRequired('Method'); % Required input 'Method'
Analyse_Var.addRequired('Geometry2D'); % Required input 'Geometry2D'
Analyse_Var.addRequired('Physics2D'); % Required input 'Physics2D'
Analyse_Var.addRequired('Rabi'); % Required input 'Rabi' with default value 1
Analyse_Var.addRequired('delta'); % Required input 'delta' with default value 0
%Analyse_Var.addOptional('InitialData_Choice',1); % Optional input 'InitialData_Choice' with default value 1
Analyse_Var.addOptional('X0', 0); % Optional input 'X0' with default value 0
Analyse_Var.addOptional('Y0', 0); % Optional input 'Y0' with default value 0
Analyse_Var.addOptional('gamma_x', 1); % Optional input 'gamma_x' with default value 1
Analyse_Var.addOptional('gamma_y', 1); % Optional input 'gamma_y' with default value 1

%% Parsing inputs and storing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
Method = Analyse_Var.Results.Method; %Storing the 'Method' input
Geometry2D = Analyse_Var.Results.Geometry2D; %Storing the 'Geometry2D' input
Physics2D = Analyse_Var.Results.Physics2D; %Storing the 'Physics2D' input
%InitialData_Choice = Analyse_Var.Results.InitialData_Choice; %Storing the 'InitialData_Choice' input
X0 = Analyse_Var.Results.X0; %Storing the 'X0' input
% IF X0 is a scalar
if (isscalar(X0))
   X0 = X0*ones(1,Method.Ncomponents); % Storing 'X0' as a vector
end
Y0 = Analyse_Var.Results.Y0; %Storing the 'Y0' input
% IF Y0 is a scalar
if (isscalar(Y0))
   Y0 = Y0*ones(1,Method.Ncomponents); % Storing 'Y0' as a vector
end
gamma_x = Analyse_Var.Results.gamma_x; %Storing the 'gamma_x' input
gamma_y = Analyse_Var.Results.gamma_y; %Storing the 'gamma_y' input

Rabi = Analyse_Var.Results.Rabi; %Storing the 'Omega' input
delta = Analyse_Var.Results.delta; %Storing the 'delta' input

%% Initialization of the initial data variable
Phi_0 = cell(1,Method.Ncomponents); % Initialization of the initial data

% Compute wavefunctions with correct density weights
theta = pi/2-atan(delta/Rabi); % Compute theta from cotangent relation
rho_up = sin(theta/2)^2;        % Probability weight for spin-up
rho_down = cos(theta/2)^2;      % Probability weight for spin-down

A = (gamma_x * gamma_y)^(1/4) / sqrt(pi); % Normalization factor
Phi_0{1} = A.*exp(-(gamma_x*(Geometry2D.X - X0(1)).^2 + gamma_y*(Geometry2D.Y - Y0(1)).^2)/2) + 1i * zeros(size(Geometry2D.X));   % Spin-up component
Phi_0{2} = A.*exp(-(gamma_x*(Geometry2D.X - X0(2)).^2 + gamma_y*(Geometry2D.Y - Y0(2)).^2)/2) + 1i * zeros(size(Geometry2D.X)); % Spin-down component

% I could also use the custom function for gaussian profiles
% Phi_0{1} = GaussianInitialData2d(Geometry2D, Physics2D, gamma_x, gamma_y, X0(1), Y0(1));
% Phi_0{2} = GaussianInitialData2d(Geometry2D, Physics2D, gamma_x, gamma_y, X0(2), Y0(2));

% Normalize each component separately
Phi_0{1} = sqrt(rho_up) * Phi_0{1} / L2_norm2d(Phi_0{1}, Geometry2D);
Phi_0{2} = sqrt(rho_down) * Phi_0{2} / L2_norm2d(Phi_0{2}, Geometry2D); 
   
end