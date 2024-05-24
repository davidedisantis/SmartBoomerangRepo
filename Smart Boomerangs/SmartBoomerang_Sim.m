%% Data
clear
% close all
clc

tic
%--------------------------Simulation parameters--------------------------%
t = 10;                     % Simulation time [s]
dt = 2e-4;                    % Fixed-step timestep
f_con = 1000;                 % Control frequency [Hz]                
dt_con = 1/f_con;             % Control timestep [s]
minstep = 1e-5;               % Min timestep for variable-step simulation [s]
maxstep = 1.66e-4;            % Max timestep for variable-step simulation [s]

%-----------------------Boomerang design parameters-----------------------%
l_blade = 25e-2;       % Length of one blade [m]
nw = 4;                 % Number of wings [ ]
R = l_blade;            % Radius of rotation [m]
S = pi*R^2;             % Disk area [m^2]
m = 200e-3;          % Boomerang mass [kg]
c = 5e-2;               % Mean chord [m]
S_b = l_blade*c*nw;

d_ac = c/4; % Position of aerodynamic center in body coordinates
LAMBDAj = 0;  % Wing sweep angle [deg]
LAMBDAj = deg2rad(LAMBDAj);
gamma = 90;    % Folding angle [deg]
gamma = deg2rad(gamma);
betaj = 4;      % Wing coning angle 
betaj = deg2rad(betaj);
thetaj = 0;     % Wing pitch angle at root
x_ac = zeros(3,nw);
for i = 1:nw
    Rac = [cos(gamma*i), sin(gamma*i), 0; -sin(gamma*i), cos(gamma*i), 0;0, 0, 1];
    x_ac(:,i) = Rac*[0; -d_ac; 0];
end
%-------------------------Launch conditions-------------------------------%
ThAng = 45;                                 % Throw angle from East direction
U0mod = 25;
U0 = [U0mod*cosd(ThAng); U0mod*sind(ThAng); 0];   % Initial throw speed in inertial frame [m/s]
omega0 = [0; 0; 10*pi*2];                   % Initial angular speed [rad/s]
PHI0 = 75;                                  % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
THETA0 = 0;                                 % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);  
PSI0 = 220;                                 % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
R_pos0 = [0; 0; 1.5];                       % Initial position in inertial frame
%-----------------Moments of inertia of a single blade--------------------%
I_xi = (m/nw)/12*(l_blade^2 + (0.1*c)^2) + m/nw*(l_blade/2)^2; 
I_eta = (m/nw)/12*(c^2);
I_zeta = I_xi+I_eta;%(m_b/nw)/12*(l_blade^2 + c^2) + m_b/nw*(l_blade/2)^2;
I_xieta = 0;

Jj = [
    I_xi I_xieta 0;
    I_xieta I_eta 0;
    0 0 I_zeta
    ];
%---------------------------Planetary parameters--------------------------%
g = 9.81;                % Gravity acceleration [m/s^2]
rho = 1.22;              % Air density at sea level [kg/m^3]
v_wind = [0; 1.5; 0];      % Wind speed in inertial reference frame [m/s]
% v_wind = [-1.5; 0; 0];

%------------Parameters for integration of aerodynamic forces-------------%
n = 50; % Number of intervals
n = n + 1; 
l_integrate = linspace(0,l_blade,n);
%---------------------------Rotation matrices-----------------------------%
Tj = zeros(3,3, nw);
invTj = zeros(3,3,nw);

for i = 0:(nw-1)
    Rj1 = [1, 0, 0; 0, cos(betaj), sin(betaj); 0, -sin(betaj), cos(betaj)];
    Rj2 = [cos(thetaj), 0, -sin(thetaj); 0, 1, 0; sin(thetaj), 0, cos(thetaj)];
    Rj3 = [cos(LAMBDAj-pi/2 + gamma*i), sin(LAMBDAj-pi/2 + gamma*i), 0; ...
        -sin(LAMBDAj-pi/2 + gamma*i), cos(LAMBDAj-pi/2+ gamma*i), 0;0, 0, 1];
    Tj(:,:,i+1) = Rj2*Rj1*Rj3;
    invTj(:,:,i+1) = transpose(Tj(:,:,i+1));
end
RI01 = [1, 0, 0; 0, cos(PHI0), sin(PHI0); 0, -sin(PHI0), cos(PHI0)];
RI02 = [cos(THETA0), 0, -sin(THETA0); 0, 1, 0; sin(THETA0), 0, cos(THETA0)];
RI03 = [cos(PSI0), sin(PSI0), 0; -sin(PSI0), cos(PSI0), 0; 0, 0, 1];
TI0 = RI01*RI02*RI03;
invTI0 = transpose(TI0);

%------------Computation of moment of inertia wrt body frame--------------%
Ji = zeros(3,3,nw);
J = zeros(3);

for i = 1:nw
    Ji(:,:,i) = invTj(:,:,i)*Jj*Tj(:,:,i);
    J = J + (Ji(:,:,i));
    % for j = 1:3
    %     for k = 1:3
    %         J(j,k) = J(j,k) + m_b/nw*(norm(x_ac(:,nw))^2 - x_ac(j,nw)*x_ac(k,nw));
    %     end
    % end
end

% J(2,2) = J(2,2) + m_b*x_ac^2;
% J(3,3) = J(3,3) + m_b*x_ac^2;

% J = diag(diag(J));

invJ = inv(J);

%-------------------Starting conditions in body frame---------------------%
Eul0 = [0; 0; 0];
Quat0 = [1; 0; 0; 0];
u0 = TI0*U0;
r0 = TI0*R_pos0;

%------------------------Aerodynamic coefficients-------------------------%
% CL data
M1 = readmatrix("Cl.csv");
x_CL = M1(:,3);
y_CL = 1.35*M1(:,2);
% y_CL = M1(:,2);
CLdata = [x_CL, y_CL];

% CD data
M2 = readmatrix("Cd.csv");
x_CD = M2(:,3);
y_CD = 1.35*M2(:,2);
% y_CD = M2(:,2);
CDdata = [x_CD, y_CD];

% CM data
M3 = readmatrix("Cm.csv");
x_CM = M3(:,3);
y_CM = M3(:,2);
CMdata = [x_CM, y_CM];


%---------------------Desired trajectory: returning path------------------% 
load('Returning_Traj_4wings.mat');
x_des = R_inertial1(:,1);
y_des = R_inertial1(:,2);

X_des = [x_des, y_des];

choice = input(['Enter 1 to choose a wing-pitch+attitude control method\n' ...
                'Enter 2 to choose a wing-pitch control method\n' ...
                'Enter 3 to choose a dihedral angle control method\n' ...
                'Enter 4 to chose a blade-length control method\n']);

switch choice
    case 1
    %-------------------------Pitch control parameters------------------------%
    servo_rpm = 500;           % Rotation speed of the servo/rotation speed of the wing
    servo_rpm = deg2rad(servo_rpm);
    theta_sat_u = 40;               % Max achievable collective pitch angle [deg]
    theta_sat_l = -5;               % Min achievable collective pitch angle [deg]
    theta_sat_u = deg2rad(theta_sat_u);
    theta_sat_l = deg2rad(theta_sat_l);
    
    theta_sat_u_cy = 10;                % Max achievable cyclice pitch angle [deg]
    theta_sat_l_cy = -10;               % Min achievable cyclice pitch angle [deg]
    theta_sat_u_cy = deg2rad(theta_sat_u_cy);
    theta_sat_l_cy = deg2rad(theta_sat_l_cy);
    
    
    %-------------------------Pitch PID parameters----------------------------%
    Kp = .6; 
    Ki = 0;
    Kd = 1;
    
    krchi = 4; % k that allows to control roll to cancel error in CHI
    krz = -1;   % k that allows to control roll to cancel error in Z
    kpchi = 10; % k that allows to control wing-pitch to cancel error in CHI
    kpz = 5;
    %---------------Cyclic pitch - attitude control PID parameters------------%
    Kp_roll = 7;
    Ki_roll = 0;
    Kd_roll = 0;
    
    Kp_pitch = 8*0;
    Ki_pitch = 0;
    Kd_pitch = 2*0;
    sim = sim("Boomerang_pitch_att_simulink.slx");
    case 2
        %-------------------------Pitch control parameters------------------------%
        servo_rpm = 500;           % Rotation speed of the servo/rotation speed of the wing
        servo_rpm = deg2rad(servo_rpm);
        theta_sat_u = 40;               % Max achievable collective pitch angle [deg]
        theta_sat_l = -5;               % Min achievable collective pitch angle [deg]
        theta_sat_u = deg2rad(theta_sat_u);
        theta_sat_l = deg2rad(theta_sat_l);
        
        theta_sat_u_cy = 10;                % Max achievable cyclice pitch angle [deg]
        theta_sat_l_cy = -10;               % Min achievable cyclice pitch angle [deg]
        theta_sat_u_cy = deg2rad(theta_sat_u_cy);
        theta_sat_l_cy = deg2rad(theta_sat_l_cy);
        
        
        %-------------------------Pitch PID parameters----------------------------%
        Kp = .6; 
        Ki = 0;
        Kd = 1;
        sim = sim("Boomerang_pitch_simulink");

    case 3
        %-----------------------Dihedral control parameters-----------------------%
        servo_rpm = 500;           % Rotation speed of the servo/rotation speed of the wing
        servo_rpm = deg2rad(servo_rpm);
        max_dd = 10;              % Max delta of joint angle [deg]
        max_dd = deg2rad(max_dd);  
        min_dd = 0;              % Min delta of joint angle [deg]
        min_dd = deg2rad(min_dd);
        
        %---------------------------Dihedral PD parameters------------------------%
        Kp = 0.5;
        Ki = 0;
        Kd = 0.2;
        Kp = 0; Ki = 0; Kd = 0;
        sim = sim("Boomerang_dihedral_simulink");

    case 4
        %------------------------Blade length PID parameters----------------------%
        Kp = 0.5;
        Ki = 0;
        Kd = 1;
        
        max_l = 5e-2;               % Max positive delta of blade length
        min_l = -5e-2;              % Max negative delta of blade length
        l_rate = .2;                % Max blade length variation speed [m/s]
        sim = sim("Boomerang_lblade_simulink");
end

toc
%% Plot top view trajectory and compare with desired and uncontrolled trajectories
% u = zeros(length(sim.tout),3);
% normU = zeros(length(sim.tout),1);

%-Load the correct file for uncontrolled trajectory under effects of wind-%
if v_wind(1) == -1.5 && v_wind(2) == 0 && v_wind(3) == 0
    load('Wind_wx_4wings.mat')
else if v_wind(1) == 1.5 && v_wind(2) == 0 && v_wind(3) == 0
    load('Wind_wxx_4wings.mat')
else if v_wind(1) == -3 && v_wind(2) == 0 && v_wind(3) == 0
    load('Wind_2wx_4wings.mat')
else if v_wind(2) == 1.5 && v_wind(1) == 0 && v_wind(3) == 0
    load('Wind_wy_4wings.mat')
else if v_wind(2) == -1.5 && v_wind(1) == 0 && v_wind(3) == 0
    load('Wind_wyy_4wings.mat')
else
    R_inertial_w = [0; 0]; % Do not plot uncontrolled trajectory if it was not previously computed
end
end
end
end
end
x_w = R_inertial_w(1,:);
y_w = R_inertial_w(2,:);

R_inertial = sim.R_inertial(:,:);
fprintf('\nThe final distance from the thrower is %f m\n',norm(R_inertial(1:2,end)))

% Plot top-view trajectory 
hold on
figure(1)
plot(R_inertial(1,1),R_inertial(2,1), 'go', ...
    R_inertial(1,:), R_inertial(2,:), 'k-', ...
    R_inertial(1,end), R_inertial(2,end), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
hold on
plot(x_des,y_des, 'b--');
hold on
plot(x_w, y_w, 'r--');
legend('', 'Controlled trajectory', '', 'Desired trajectory','Uncontrolled trajectory', 'Location','southeast','FontSize',11)
axis equal
err_vec = zeros(length(R_inertial),1);

% for i = 1:length(R_inertial)
%     err_vec(i) = min(vecnorm([R_inertial(1,i) R_inertial(2,i)] - X_des,2,2));
%     % err_sign = sign(norm([R_inertial(1,i) R_inertial(2,i)]) - norm([X_des(i,1), X_des(i,2)]));
%     % err_vec(i) = err_vec(i)*err_sign;
% end
% figure(2)
% plot(sim.tout, err_vec)
%% Plot results 

% phi = sim.EulAng(:,1);
% theta = sim.EulAng(:,2);
% psi = sim.EulAng(:,3);
R_inertial = sim.R_inertial(:,:);
U = sim.U(:,:)';
normU = vecnorm(U');

% Plots
figure(1)
plot(sim.tout, R_inertial(3,:), 'k-')
xlabel('Time [s]')
ylabel('Altitude [m]')

figure(2)
plot(R_inertial(1,1),R_inertial(2,1), 'go', ...
    R_inertial(1,:), R_inertial(2,:), 'k-', ...
    R_inertial(1,end), R_inertial(2,end), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
axis equal
% xlim([-10 10])
% ylim([-5 18])

figure(3)
plot(R_inertial(1,1),R_inertial(3,1), 'go', ...
    R_inertial(1,:), R_inertial(3,:), 'k-', ...
    R_inertial(1,end), R_inertial(3,end), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
axis equal
% xlim([-10 10])
% ylim([-2 8])

figure(4)
plot(R_inertial(2,1), R_inertial(3,1), 'go', ...
    R_inertial(2,:), R_inertial(3,:), 'k-', ...
    R_inertial(2,end), R_inertial(3,end), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')
axis equal
% xlim([-5 18])
% ylim([-2 8])

figure(5)
plot(sim.tout, normU, 'k-')
xlabel('Time [s]')
ylabel('Speed U [m/s]')
% xlim([0 5])
% ylim([0 25])

% figure(6)
% plot(sim.tout, sim.omega(:,3),'k')

figure(6)
plot3(R_inertial(1,1), R_inertial(2,1), R_inertial(3,1), 'go', ...
    R_inertial(1,:), R_inertial(2,:), R_inertial(3,:), 'k-', ...
    R_inertial(1,end), R_inertial(2,end), R_inertial(3,end), 'ro')
title('3D Trajectory')
grid on

%% Trajectory and attitude - Euler Angles
phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
R_inertial = sim.R_inertial(:,:)';
% % Define the filename for the video
% videoFilename = 'trajectory_video.avi';

% % Create a VideoWriter object
% video = VideoWriter(videoFilename);
% 
% % Set the frame rate (frames per second)
% numFrames = length(phi(1:50:end));
% frameRate = numFrames/(sim.tout(end)); % Adjust as needed
% video.FrameRate = frameRate;
% 
% % Open the VideoWriter
% open(video);

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
% R_inertial = zeros(size(R_inertial));
for i = 1:50:length(sim.tout)
    R1 = [
    1, 0, 0;
    0, cos(phi(i)), sin(phi(i));
    0, -sin(phi(i)), cos(phi(i))
    ];
    R2 = [
    cos(theta(i)), 0, -sin(theta(i));
    0, 1, 0;
    sin(theta(i)), 0, cos(theta(i))
    ];
    R3 = [
    cos(psi(i)), sin(psi(i)), 0;
    -sin(psi(i)), cos(psi(i)), 0;
    0, 0, 1
    ];
    T0 = R1*R2*R3;
    invT0 = transpose(T0);
    attx = invTI0*invT0*x_v;
    atty = invTI0*invT0*y_v;
    attz = invTI0*invT0*z_v;

    figure(1);
    plot3(R_inertial(:,1),R_inertial(:,2),R_inertial(:,3),'k-')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attx(1), attx(2), attx(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), atty(1), atty(2), atty(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attz(1), attz(2), attz(3), 'r')    
    xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
    hold off
    legend('Trajectory', 'x_b', 'y_b', 'z_b','interpreter', 'TeX')
    title('Time elapsed: ', num2str(sim.tout(i)))
    grid on

    xlim([min(R_inertial(:,1))-1 max(R_inertial(:,1))+1])
    ylim([min(R_inertial(:,2))-1 max(R_inertial(:,2))+1])
    zlim([min(R_inertial(:,3))-1 max(R_inertial(:,3))+1])
    %  % Capture the current frame
    % frame = getframe(gcf);
    % 
    % % Write the frame to the video
    % writeVideo(video, frame);
end
% % Close the VideoWriter
% close(video);
%% Plot and compare energies with undisturbed returning trajectory

figure(1)
plot(sim.P_E, 'k-')
hold on
plot(P_Enc + 0*m*g, 'k--') % Adding 3.5*m*g to shift the curve as starting at 5.3 m 
title('Potential Energy over time')
xlabel('Time [s]')
ylabel('Energy [J]')
legend('E_P control', 'E_P no control', 'FontSize', 11)

figure(2)
plot(sim.K_E, 'k-')
hold on
plot(K_Enc, 'k--')
title('Kinetic Energy over time')
xlabel('Time [s]')
ylabel('Energy [J]')
legend('E_K control', 'E_K no control', 'FontSize', 11)

figure(3)
plot(sim.R_E, 'k-')
hold on
plot(R_Enc, 'k--')
title('Rotational Energy over time')
xlabel('Time [s]')
ylabel('Energy [J]')
legend('E_R control', 'E_R no control', 'FontSize', 11)

%% Trajectory and attitude of the non-spinning disk w/ Tn

phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
PHI = zeros(size(phi));
THETA = zeros(size(phi));
PSI = zeros(size(phi));
lambda = sim.lambda;
R_inertial = sim.R_inertial(:,:)';
% % Define the filename for the video
% videoFilename = 'trajectory_video_lambda.avi';

% Create a VideoWriter object
% video = VideoWriter(videoFilename);

% % Set the frame rate (frames per second)
% numFrames = length(phi(1:50:end));
% t_video = 20;                 % Desired video duration [s]
% frameRate = numFrames/t_video; % Adjust as needed
% video.FrameRate = frameRate;

% Open the VideoWriter
% open(video);

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
% R_inertial = zeros(size(R_inertial));
for i = 1:length(sim.tout)
    R1 = [1, 0, 0; 0, cos(phi(i)), sin(phi(i)); 0, -sin(phi(i)), cos(phi(i))];
    R2 = [cos(theta(i)), 0, -sin(theta(i)); 0, 1, 0; sin(theta(i)), 0, cos(theta(i))];
    R3 = [cos(psi(i)), sin(psi(i)), 0; -sin(psi(i)), cos(psi(i)), 0; 0, 0, 1];
    T0 = R1*R2*R3;
    Tn = [cos(lambda(i)) -sin(lambda(i)) 0; sin(lambda(i)) cos(lambda(i)) 0; 0 0 1];
    TI = Tn*T0*TI0;
    THETA(i) = -asin(TI(1,3));
    PSI(i) = acos( TI(1,1)/cos(THETA(i)) );
    PHI(i) = acos( TI(3,3)/cos(THETA(i)) );
end

figure(2)
plot(sim.tout, rad2deg(PHI),'LineWidth', 1)
title('Roll angle')
figure(3)
plot(sim.tout, rad2deg(THETA),'LineWidth', 1)
title('Pitch angle')
% figure(4)
% plot(sim.tout, rad2deg(PSI),'LineWidth', 1)

for i = 1:50:length(sim.tout)
    R1 = [1, 0, 0; 0, cos(phi(i)), sin(phi(i)); 0, -sin(phi(i)), cos(phi(i))];
    R2 = [cos(theta(i)), 0, -sin(theta(i)); 0, 1, 0; sin(theta(i)), 0, cos(theta(i))];
    R3 = [cos(psi(i)), sin(psi(i)), 0; -sin(psi(i)), cos(psi(i)), 0; 0, 0, 1];
    T0 = R1*R2*R3;
    invT0 = transpose(T0);
    attx = invTI0*invT0*x_v;
    atty = invTI0*invT0*y_v;
    attz = invTI0*invT0*z_v;
    
    Tn = [cos(lambda(i)) -sin(lambda(i)) 0; sin(lambda(i)) cos(lambda(i)) 0; 0 0 1];
    invTn = transpose(Tn);
    attxn = invTI0*invT0*invTn*x_v;
    attyn = invTI0*invT0*invTn*y_v;
    attzn = invTI0*invT0*invTn*z_v;
    TI = Tn*T0*TI0;
    THETA(i) = -asin(TI(1,3));
    PSI(i) = acos( TI(1,1)/cos(THETA(i)) );
    PHI(i) = acos( TI(3,3)/cos(THETA(i)) );

    figure(1);
    plot3(R_inertial(:,1),R_inertial(:,2),R_inertial(:,3),'k-')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attx(1), attx(2), attx(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), atty(1), atty(2), atty(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attz(1), attz(2), attz(3), 'r')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attxn(1), attxn(2), attxn(3), 'g')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attyn(1), attyn(2), attyn(3), 'g')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attzn(1), attzn(2), attzn(3), 'm')    
    xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
    hold off
    legend('Trajectory', 'x_b', 'y_b', 'z_b','interpreter', 'TeX')
    title('Time elapsed: ', num2str(sim.tout(i)))
    grid on

    xlim([min(R_inertial(:,1))-1 max(R_inertial(:,1))+1])
    ylim([min(R_inertial(:,2))-1 max(R_inertial(:,2))+1])
    zlim([min(R_inertial(:,3))-1 max(R_inertial(:,3))+1])
     % Capture the current frame
    % frame = getframe(gcf);
    % 
    % % Write the frame to the video
    % writeVideo(video, frame);
end
% % Close the VideoWriter
% close(video);
%% Plot roll angle and desired roll angle
close all
clc
phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
THETA = zeros(size(phi));
PHI = zeros(size(phi));
lambda = sim.lambda;
load('PHI_natural')
for i = 1:length(sim.tout)
    R1 = [1, 0, 0; 0, cos(phi(i)), sin(phi(i)); 0, -sin(phi(i)), cos(phi(i))];
    R2 = [cos(theta(i)), 0, -sin(theta(i)); 0, 1, 0; sin(theta(i)), 0, cos(theta(i))];
    R3 = [cos(psi(i)), sin(psi(i)), 0; -sin(psi(i)), cos(psi(i)), 0; 0, 0, 1];
    T0 = R1*R2*R3;
    Tn = [cos(lambda(i)) -sin(lambda(i)) 0; sin(lambda(i)) cos(lambda(i)) 0; 0 0 1];
    TI = Tn*T0*TI0;
    THETA(i) = -asin(TI(1,3));
    PSI(i) = acos( TI(1,1)/cos(THETA(i)) );
    PHI(i) = acos( TI(3,3)/cos(THETA(i)) );
end
% figure(2)
plot(sim.tout, rad2deg(PHI), sim.tout, rad2deg(PHI_des)*ones(size(sim.tout)),'LineWidth', 1)
hold on
plot(t1, rad2deg(PHI1), 'k--')
legend('Roll angle', 'Desired roll angle', 'Uncontrolled roll angle')
% figure(3)
% plot(sim.tout, rad2deg(THETA),'LineWidth', 1)
% figure(4)
% plot(sim.tout, rad2deg(PSI),'LineWidth', 1)

% Averages
[~,maxind] = findpeaks(lambda); % Find the local maxima of lambda 
maxind = [1;maxind];
mean_t = zeros(length(maxind)-1,1);
mean_THETA = zeros(length(maxind)-1,1);

for i = 1:(length(maxind)-1)
    mean_t(i) = mean(sim.tout(maxind(i):maxind(i+1)));
    mean_THETA(i) = mean(PHI(maxind(i):maxind(i+1)));
end

figure(2)
plot(sim.tout,rad2deg(PHI_des)*ones(size(sim.tout)), mean_t, rad2deg(mean_THETA), 'LineWidth',1.5);
hold on
plot(mean_t1, rad2deg(mean_PHI1), 'k--','LineWidth', 0.5)
xlim([0, sim.tout(end)])
legend('Desired angle', 'Avg Actual angle', 'Avg Uncontrolled angle - natural trajectory', 'FontSize', 11)
xlabel('Time [s]')
ylabel('Roll angle \Phi [deg]')

s = abs(max(abs(mean_THETA)) - PHI_des);
fprintf('\nThe maximum overshoot is %f degrees \n', rad2deg(s))

%% Plot roll angle and desired pitch angle % da rimuovere in versione finale, serve solo come prova 
close all
clc
THETA_des = 5*ones(length(sim.tout),1);

phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
THETA = zeros(size(phi));
PHI = zeros(size(phi));
lambda = sim.lambda;
load('PHI_natural')
for i = 1:length(sim.tout)
    R1 = [1, 0, 0; 0, cos(phi(i)), sin(phi(i)); 0, -sin(phi(i)), cos(phi(i))];
    R2 = [cos(theta(i)), 0, -sin(theta(i)); 0, 1, 0; sin(theta(i)), 0, cos(theta(i))];
    R3 = [cos(psi(i)), sin(psi(i)), 0; -sin(psi(i)), cos(psi(i)), 0; 0, 0, 1];
    T0 = R1*R2*R3;
    Tn = [cos(lambda(i)) -sin(lambda(i)) 0; sin(lambda(i)) cos(lambda(i)) 0; 0 0 1];
    TI = Tn*T0*TI0;
    THETA(i) = -asin(TI(1,3));
    PSI(i) = acos( TI(1,1)/cos(THETA(i)) );
    PHI(i) = acos( TI(3,3)/cos(THETA(i)) );
end
% figure(2)
plot(sim.tout, rad2deg(THETA), sim.tout, (THETA_des),'LineWidth', 1)
hold on
% plot(t1, rad2deg(PHI1), 'k--')
% legend('Roll angle', 'Desired roll angle', 'Uncontrolled roll angle')
% figure(3)
% plot(sim.tout, rad2deg(THETA),'LineWidth', 1)
% figure(4)
% plot(sim.tout, rad2deg(PSI),'LineWidth', 1)

% Averages
[~,maxind] = findpeaks(lambda); % Find the local maxima of lambda 
maxind = [1;maxind];
mean_t = zeros(length(maxind)-1,1);
mean_THETA = zeros(length(maxind)-1,1);

for i = 1:(length(maxind)-1)
    mean_t(i) = mean(sim.tout(maxind(i):maxind(i+1)));
    mean_THETA(i) = mean(THETA(maxind(i):maxind(i+1)));
end
figure(2)
plot(sim.tout,THETA_des, mean_t, rad2deg(mean_THETA), 'LineWidth',1.5);
% hold on
% plot(mean_t1, rad2deg(mean_PHI1), 'k--','LineWidth', 0.5)
xlim([0, sim.tout(end)])
% legend('Desired angle', 'Avg Actual angle', 'Avg Uncontrolled angle - natural trajectory', 'FontSize', 11)
% xlabel('Time [s]')
% ylabel('Roll angle \Phi [deg]')

% s = abs(max(abs(mean_THETA)) - PHI_des);
% fprintf('\nThe maximum overshoot is %f degrees \n', rad2deg(s))

%% Energy computation test
% EW = sim.EW; % Wind energy
% ED = sim.ED; % Aerodynamic energy
% U = sim.U(:,:);
% normU = vecnorm(U);
% R_inertial = sim.R_inertial(:,:);
% 
% plot(sim.tout, EW + ED - 1/2*m*normU'.^2 - m*g*R_inertial(3,:)')
% hold on
% plot(sim.tout, EW)

%% Plot results - cut when close to return point

phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
R_inertial = sim.R_inertial(:,:);
U = sim.U(:,:)';
normU = vecnorm(U');

[~,ind] = min((vecnorm(R_inertial(1:2,ceil(length(R_inertial)/2):end),2)));

ind = ceil(length(R_inertial)/2) + ind - 1;
R_inertial = R_inertial(:,1:ind);
t = sim.tout(1:ind);

% Plots
figure(1)
plot(t, R_inertial(3,:), 'k-')
xlabel('Time [s]')
ylabel('Altitude [m]')

figure(2)
plot(R_inertial(1,1),R_inertial(2,1), 'go', ...
    R_inertial(1,:), R_inertial(2,:), 'k-', ...
    R_inertial(1,end), R_inertial(2,end), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
axis equal
% xlim([-10 10])
% ylim([-5 18])

figure(3)
plot(R_inertial(1,1),R_inertial(3,1), 'go', ...
    R_inertial(1,:), R_inertial(3,:), 'k-', ...
    R_inertial(1,end), R_inertial(3,end), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
axis equal
% xlim([-10 10])
% ylim([-2 8])

figure(4)
plot(R_inertial(2,1), R_inertial(3,1), 'go', ...
    R_inertial(2,:), R_inertial(3,:), 'k-', ...
    R_inertial(2,end), R_inertial(3,end), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')
axis equal
% xlim([-5 18])
% ylim([-2 8])

figure(5)
plot(sim.tout, normU, 'k-')
xlabel('Time [s]')
ylabel('Speed U [m/s]')
% xlim([0 5])
% ylim([0 25])

% figure(6)
% plot(sim.tout, sim.omega(:,3),'k')

figure(6)
plot3(R_inertial(1,1), R_inertial(2,1), R_inertial(3,1), 'go', ...
    R_inertial(1,:), R_inertial(2,:), R_inertial(3,:), 'k-', ...
    R_inertial(1,end), R_inertial(2,end), R_inertial(3,end), 'ro')
title('3D Trajectory')
grid on

%% Plot control results
clear
clc
close all

prompt = ['Enter 1 to visualize controlled trajectory with disturbance of a wind blowing to the west\n' ...
    'Enter 2 to visualize controlled trajectory with disturbance of a wind blowing to the south\n' ...
    'Input: '];
answer = input(prompt);
load('Returning_Traj_4wings.mat')
x_des = R_inertial1(:,1);
y_des = R_inertial1(:,2);

if answer == 1
    load('Control_-wx.mat');
    load('Wind_wx_4wings.mat');
else 
    if answer == 2
        load('Control_-wy.mat');
        load('Wind_wyy_4wings.mat');
    else 
        printf('\nInput Error\n')
    end
end

x_w = R_inertial_w(1,:);
y_w = R_inertial_w(2,:);

% Plots
figure(1)
plot(t, R_inertial(3,:), 'k-')
xlabel('Time [s]')
ylabel('Altitude [m]')

figure(2)
plot(R_inertial(1,1),R_inertial(2,1), 'go', ...
    R_inertial(1,:), R_inertial(2,:), 'k-', ...
    R_inertial(1,end), R_inertial(2,end), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
hold on
plot(x_des,y_des, 'b--');
hold on
plot(x_w, y_w, 'r--');
legend('', 'Controlled trajectory', '', 'Desired trajectory','Uncontrolled trajectory', 'Location','southeast','FontSize',11)
axis equal
% xlim([-10 10])
% ylim([-5 18])

figure(3)
plot(R_inertial(1,1),R_inertial(3,1), 'go', ...
    R_inertial(1,:), R_inertial(3,:), 'k-', ...
    R_inertial(1,end), R_inertial(3,end), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
axis equal
% xlim([-10 10])
% ylim([-2 8])

figure(4)
plot(R_inertial(2,1), R_inertial(3,1), 'go', ...
    R_inertial(2,:), R_inertial(3,:), 'k-', ...
    R_inertial(2,end), R_inertial(3,end), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')
axis equal
% xlim([-5 18])
% ylim([-2 8])

figure(5)
plot(t, normU, 'k-')
xlabel('Time [s]')
ylabel('Speed U [m/s]')
% xlim([0 5])
% ylim([0 25])

% figure(6)
% plot(t, sim.omega(:,3),'k')

figure(6)
plot3(R_inertial(1,1), R_inertial(2,1), R_inertial(3,1), 'go', ...
    R_inertial(1,:), R_inertial(2,:), R_inertial(3,:), 'k-', ...
    R_inertial(1,end), R_inertial(2,end), R_inertial(3,end), 'ro')
title('3D Trajectory')
grid on