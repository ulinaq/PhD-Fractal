function prm = Tree_parameter(alpha)
%% Artery Parameters
% blood and flow parameters (all unit in cgs)
prm.Rho = 1;
prm.mu = 1;
prm.Pin = 50;              % capillary pressure
prm.Qin = 1;               % Flowrate in input node

% Geometry parameters
prm.Ro = alpha;        % Initial radius 
prm.LtoR = 9;        % Ratio between length and radius
prm.Mexp = 2.7;       % Bifurcation exponent
prm.Rmin = 0.8*prm.Ro;      % Minimum radius
prm.N = 2;            % number of bifurcation
gamma = 1;          % Radius ratio between both childs vessel (degree of symmetry)

% Ratio between parent and 2 childs vessel
delta = (1+gamma^(prm.Mexp/2))^(-1/prm.Mexp); 
beta = delta*sqrt(gamma);
prm.phi(1) = acos((1+delta^4-beta^4)/(2*delta^2));    % angle between first and parent
prm.phi(2) = acos((1-delta^4+beta^4)/(2*beta^2));     % angle between second and parent
prm.ratio = [delta, beta];   % Radius ratio between 2 childs to the parent vessel

%% Vein Parameters
% Geometry parameters
prm.Rvo = 1.5*prm.Ro;       % Initial radius 
prm.Mexp = 2.8;             % Bifurcation exponent
prm.Rvmin = 1.5*prm.Rmin;   % Minimum radius
prm.N = 2;                  % number of bifurcation
gamma = 1;                  % Radius ratio between both childs vessel (degree of symmetry)

% Ratio between parent and 2 childs vessel
delta = (1+gamma^(prm.Mexp/2))^(-1/prm.Mexp); 
beta = delta*sqrt(gamma);
prm.vphi(1) = acos((1+delta^4-beta^4)/(2*delta^2));    % angle between first and parent
prm.vphi(2) = acos((1-delta^4+beta^4)/(2*beta^2));     % angle between second and parent
prm.ratio = [delta, beta];   % Radius ratio between 2 childs to the parent vessel
end