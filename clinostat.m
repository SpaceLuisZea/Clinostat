%{
FILENAME:  clinostat

DESCRIPTION: This script is used to investigate the appropriate clinostat
rotational speed with the following objectives in mind:

(I) to determine the range of angular velocities that will produce a "q"
value less than 10 microns;

(II) to determine how long it will take for a particle to exit the
deplection zone. 

This work is based on the analysis carred out in:
KESSLER, J.O., "The Internal Dynamics of Slowly Rotating Biological
Systems," ASGSB Bull. 5(2):11-21, (1992)

FUNCTIONS/SCRIPTS CALLED:
    clino_dyn

MODIFICATIONS:
    13-Oct-2018    |    Aaron Rosengren - Original
%}

% clearvars; close all; clc
tic

% -------------------------------------------------------------------------
%     CONSTANTS & PARAMETERS   
% -------------------------------------------------------------------------

%  Constants
g = 980.7;                                      % acceleration of gravity [cm/s^2]
rpm2radps = 2*pi/60;                            % conversion from revolutions per min to radians per second
sec2day = 1/86400;                              % conversion from seconds to days

%  Cellular parameters
rho = 1.08;                                     % particle density [g/cm^3]
a   = 5.599021e-5;                                % particle radius [cm]

%  Media parameters
rho_f = 1.011;                                  % fluid density [g/cm^3]
mu_f  = 0.0102;                                 % dynamic viscosity [g/cm.s]

%  Other parameters
%f     = 6*pi*mu_f*a;                            % Stokes drag force coefficient [g/s]
%m     = 3.5e-13;                            % particle mass [g]
%mstar = 5.0731E-14;                               % buoyancy-corrected mass [g]
d     = 1;                                  % depletion zone radius [cm] (0.001 for bacteria; 0.01 for fungi; 1.0 for lympho)
R     = 5;                                 % radius of vessel being clinorotated

% Read input data from file
%fileID = fopen('data_input_validation.txt','r');
fileID = fopen('data_input_small2e.txt','r');
formatSpec = '%f %f %f %f %f %f';
sizeA = [6 Inf];
data_input = fscanf(fileID,formatSpec,sizeA)';
fclose(fileID);

nrows = size(data_input,1);
ncolumns = size(data_input,2);

data_output = zeros(nrows,11);

% Loop through ICs
for z = 1:nrows

     f = data_input(z,1);
     m = data_input(z,2);
     mstar = data_input(z,3);
     x_initial = data_input(z,4);
     y_initial = data_input(z,5);
     rpm = data_input(z,6);

    omega = rpm*rpm2radps;
    state0 = [ x_initial ; y_initial ; -y_initial*omega ; -x_initial*omega ];
    r_initial = sqrt(state0(1)^2 + state0(2)^2);

    % %  KESSLER'S VALUES
    % omega = 1*rpm2radps;
    % a = 0.0003;
    % rho_f = rho-0.1;
    % R = 0.1;
    % d = 1E-02;

    % -------------------------------------------------------------------------
    %     ANALYTICAL CALCULATED PARAMETERS   
    % -------------------------------------------------------------------------

    %  Velocity of steady sedimentation [cm/s]
    v_sed = 2/9*g*(rho - rho_f)*a^2/mu_f;           % KESSLER Eqn.(1)

    %disp(' ')
    %disp('Note that the range for the sedimentation velocity for unattached')
    %disp('intracellular organelles is 1E-07 to 1E-02 cm/s.');
    %disp(' ')
    %disp(['For the parameters considered herein, we find v_sed = ', num2str(v_sed), ' cm/s.']);

    %  Exponential drift terms and predictions from analytical solution to 
    %  EOM presented by KESSLER 
    tau = g/v_sed/omega^2;                          % centrifugal time, KESSLER Eqn.(7)
    %tstar = tau*d/R*sec2day;                        % approximate time before reaching cell boundary [days]
    %tstar_dz = tau*d/sqrt((state0(1))^2+(state0(2))^2)*sec2day; % approximate time before reaching cell boundary, changing R to sqrt(x?2 + y?2) [days]
    %tstar_ev = tau*(R-sqrt((state0(1))^2+(state0(2))^2))/sqrt((state0(1))^2+(state0(2))^2)*sec2day;

    %disp(' ')
    %disp('According to the analytical approximation made by KESSLER,')
    %disp(['we expect it to take ', num2str(tau*sec2day), ' days to leave the depletion zone (tau)']);
    %disp(['we expect it to take ', num2str(tstar_dz), ' days to leave the depletion zone (tstar_dz).']);
    %disp(['we expect it to take ', num2str(tstar_ev), ' days to reach the vessel edge (tstar_ev).']);

    %  "q" radius [cm]
    q = v_sed/omega;                                % KESSLER Eqn.(13)

    %disp(' ')
    %disp('According to the analytical approximation made by KESSLER,')
    %disp(['the radius of the particle motion is predicted to be q = ', num2str(q), ' cm.']);
    %disp(' ')

    % -------------------------------------------------------------------------
    %     NUMERICALLY CALCULATED PARAMETERS   
    % -------------------------------------------------------------------------

    %  Group all constants needed in dynamics into structure
    const.g     = g;
    const.omega = omega;
    const.f     = f;
    const.m     = m;
    const.mstar = mstar;

    %  Set integrator options
    tol = 1e-9;
    options = odeset('RelTol', tol, 'AbsTol', tol );

    %  Specify timestep and integration time
    Ndays = 2;
    % Ndays = ceil(tstar);                                % could set the simulation time to estimated exit time (uncomment) 
    Nt    = 1000;
    tvec  = linspace(0, Ndays/sec2day, Nt);

    %  State propagation                 
    [time, statet] = ode15s(@clino_dyn, tvec, state0, options, const);

    radius = zeros(size(time));
    for k = 1:length(time)

        radius(k) = norm(statet(k,1:2));

    end

    tdays = (time - time(1))*sec2day;   

    x_final = statet(Nt,1);  % Final X position at end of Ndays
    y_final = statet(Nt,2);  % Final Y position at end of Ndays
    r_final = sqrt(x_final^2 + y_final^2);  % Final radius at end of Ndays

    delta_R = r_final - r_initial;   % Radial displacement during Ndays

    v_cent = delta_R / (Ndays*24*60*60);   % Centripetal speed cm/s
    %v_cent = omega * 2 * pi() * [0.5 1 2 3 4]; % Centripetal speed cm/s
    
    t_dz = d ./ v_cent;    % Time for cell to exit the depletion zone
    %t_R = (R - r_initial) ./ v_cent;     % Time for cell to to reach edge of vessel
    t_R = (R - [0.5 1 2 3 4]) ./ v_cent;     % Time for cell to to reach edge of vessel

    data_output(z,1) = q;
    data_output(z,2) = t_dz(1);
    %data_output(z,3) = t_dz(2);
    %data_output(z,4) = t_dz(3);
    %data_output(z,5) = t_dz(4);
    %data_output(z,6) = t_dz(5);
    data_output(z,7) = t_R(1);
    data_output(z,8) = t_R(2);
    data_output(z,9) = t_R(3);
    data_output(z,10) = t_R(4);
    data_output(z,11) = t_R(5);

    % hFig = figure(1); 
    %     hold all
    %     scatter( statet(:,1) , statet(:,2) , [] , tdays );
    %     colorbar
    %     hXLabel = xlabel( [ '$X', '$' ] , ...
    %                         'Interpreter', 'LaTex' );
    %     hYLabel = ylabel( [ '$Y', '$' ] , ...
    %                         'Interpreter', 'LaTex' , 'Rotation', 90 );      
    %     set( gca , 'FontName' , 'Helvetica' );
    %     set( hXLabel , 'FontSize' , 12 )
    %     set( hYLabel , 'FontSize' , 12 )
    %     set( hFig , 'units' , 'inches' , ... 
    %          'NumberTitle' , 'off' , 'Name' , 'sedimentation' );
    %     set( hFig , 'position' , [0,6.75,6,4.75] ); 

    %disp(' ')
end    
    
toc

fid = fopen('data_output.txt','wt');
for ii = 1:size(data_output,1)
    fprintf(fid,'%g\t',data_output(ii,:));
    fprintf(fid,'\n');
end
fclose(fid)
