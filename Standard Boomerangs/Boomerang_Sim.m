%% Data
clear
% close all
clc

tic
%--------------------------Simulation parameters--------------------------%
t = 2.85;                     % Simulation time [s]
dt = 2e-4;                 % Fixed-step timestep
minstep = 1e-5;            % Min timestep for variable-step simulation [s]
maxstep = 1.66e-4;         % Max timestep for variable-step simulation [s]

%-----------------------Boomerang design parameters-----------------------%
l_blade = 25e-2;       % Length of one blade [m]
nw = 4;                 % Number of wings [ ]
R = l_blade;            % Radius of rotation [m]
S = pi*R^2;             % Disk area [m^2]
m_b = 200e-3;          % Boomerang mass [kg]
c = 5e-2;               % Mean chord [m]
% m_add = zeros(3);
% m_add(3,3) = 8/3*rho*R^3;
% m = m_b*eye(3) + m_add;
m = m_b;
invm = inv(m);

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

%---------------------------Planetary parameters--------------------------%
choice = input(['Enter 1 to simulate a trajectory on Earth at sea level\n' ...
                'Enter 2 to simulate a trajectory on Titan\n' ...
                'Enter 3 to simulate a trajectory on Venus at 52 km\n' ...
                'Enter 4 to simulate a trajectory on Venus at 60 km\n']);
switch choice
    case 1
        g = 9.81;                     % Gravity acceleration on Earth[m/s^2]
        rho = 1.22;                   % Air density at sea level on Earth [kg/m^3]
    case 2
        g = 1.35;                     % Gravity acceleration on Titan [m/s^2]
        rho = 4.4*1.22;               % Air density at sea level on Titan [kg/m^3]
    case 3
        g = 8.87;                     % Gravity acceleration on Venus [m/s^2] 
        rho = 0.4854;                 % Air density at 60 km on Venus [kg/m^3]
    case 4
        g = 8.87;                     % Gravity acceleration on Venus [m/s^2] 
        rho = 1.3288;                 % Air density at 52 km on Venus [kg/m^3]
    otherwise 
        disp('Invalid input')
        check = 0;
end

v_wind = [0; 0; 0];      % Wind speed in inertial reference frame [m/s]
%-----------------Moments of inertia of a single blade--------------------%
I_xi = (m_b/nw)/12*(l_blade^2 + (0.1*c)^2) + m_b/nw*(l_blade/2)^2; 
I_eta = (m_b/nw)/12*(c^2);
I_zeta = I_xi+I_eta;%(m_b/nw)/12*(l_blade^2 + c^2) + m_b/nw*(l_blade/2)^2;
I_xieta = 0;

Jj = [
    I_xi I_xieta 0;
    I_xieta I_eta 0;
    0 0 I_zeta
    ];

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


choice2 = input(['Enter 1 to simulate a trajectory with no sensors\n' ...
                'Enter 2 to simulate a trajectory with IMU sensors\n']);
if choice2 == 1
    sim = sim("Boomerang_Simulink");
else
    %--------------------------IMU Sensor modeling----------------------------% Present data = ISM330DHCX IMU
    acc = struct('scalefactor', [], 'bias',[], 'lim', [], 'noise', []);
    gyro = struct('scalefactor', [], 'bias', [], 'gbias', [], 'lim',[], 'noise',[]);
    f = 1000;                            % Sampling frequency [Hz]
    acc.sf = [1 0.005 0.005; 0.005 1 0.005; 0.005 0.005 1]; % Scale factor and cross-couplings []
    acc.bias = [1 1 1].*60e-3*g;
    acc.lim = [-16*g  -16*g  -16*g  16*g 16*g 16*g];
    acc.noise = [80e-6*g/sqrt(f) 80e-6*g/sqrt(f) 80e-6*g/sqrt(f)];
    
    gyro.sf = [1 0.01 0.01; 0.01 1 0.01; 0.01 0.01 1];
    gbias = deg2rad(3);
    gyro.bias = [gbias gbias gbias];
    gyro.gbias = [deg2rad(0.1) deg2rad(0.1) deg2rad(0.1)];
    glim = deg2rad(4000);
    gyro.lim = [-glim  -glim  -glim  glim glim glim];
    gyro.noise = [deg2rad(7e-3)/sqrt(f) deg2rad(7e-3)/sqrt(f) deg2rad(7e-3)/sqrt(f)];
    
    noise = [acc.noise gyro.noise];
    sim = sim("Boomerang_Simulink_Nav");
end

toc
%% Plot results 
phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
R_inertial = sim.R_inertial(:,:)';
U = sim.U(:,:)';
normU = vecnorm(U');
PHI = zeros(size(phi));
THETA = zeros(size(phi));
PSI = zeros(size(phi));
lambda = sim.lambda;
rMA = zeros(3, length(sim.tout));
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
    rMA(:,i) = Tn*sim.MA(:,i);
end

% Plots
figure(1)
plot(sim.tout, rad2deg(PHI),'LineWidth', 1)
xlabel('Time [s]')
ylabel('\Phi [deg]')
hold on
plot(sim.tout, rad2deg(THETA),'LineWidth', 1)
xlabel('Time [s]')
ylabel('Angle [deg]')
legend('\Phi', '\Theta')

% figure(1)
% plot(sim.tout, R_inertial(:,3), 'k-')
% xlabel('Time [s]')
% ylabel('Altitude [m]')

figure(2)
plot(R_inertial(1,1),R_inertial(1,2), 'go', ...
    R_inertial(:,1), R_inertial(:,2), 'k-', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
axis equal
% xlim([-10 10])
% ylim([-5 18])
% xlim([-8.5 6.5])
% ylim([-3 12])

figure(3)
plot(R_inertial(1,1),R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,3), 'k-', ...
    R_inertial(end,1), R_inertial(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
axis equal
% xlim([-10 10])
% ylim([-2 8])

figure(4)
plot(R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,2), R_inertial(:,3), 'k-', ...
    R_inertial(end,2), R_inertial(end,3), 'ro')
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

figure(6)
plot(sim.tout, sim.omega(:,3),'k')
title('Spin speed')

% figure(6)
% plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
%     R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k--', ...
%     R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
% title('3D Trajectory')
% grid on

fprintf('\nFinal deviation from the thrower: %f m\n', norm([R_inertial(end,1),R_inertial(end,2)]))
fprintf('\nTotal flight time: %f s\n', sim.tout(end))
fprintf('\nMaximum height reached: %f m\n', max(R_inertial(:,3)))


%% Attitude of the non-spinning disk w/ Tn
phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
PHI = zeros(size(phi));
THETA = zeros(size(phi));
PSI = zeros(size(phi));
lambda = sim.lambda;
R_inertial = sim.R_inertial(:,:)';
% % Define the filename for the video
% videoFilename = 'trajectory_video_nospindisk.avi';

% % Create a VideoWriter object
% video = VideoWriter(videoFilename);
% 
% % Set the frame rate (frames per second)
% numFrames = length(phi(1:50:end));
% t_video = 10;                 % Desired video duration [s]
% frameRate = numFrames/t_video; % Adjust as needed
% video.FrameRate = frameRate;
% 
% % Open the VideoWriter
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
xlabel('Time [s]')
ylabel('\Phi [deg]')
hold on
plot(sim.tout, rad2deg(THETA),'LineWidth', 1)
xlabel('Time [s]')
ylabel('Angle [deg]')
legend('\Phi', '\Theta')
figure(4)
plot(sim.tout, rad2deg(PSI),'LineWidth', 1)
xlabel('Time [s]')
ylabel('\Psi [deg]')
max(rad2deg(PHI))
%% Trajectory and attitude - Euler Angles
phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
R_inertial = sim.R_inertial(:,:)';
% Define the filename for the video
videoFilename = 'trajectory_video_4w.avi';

% Create a VideoWriter object
video = VideoWriter(videoFilename);

% Set the frame rate (frames per second)
numFrames = length(phi(1:50:end));
frameRate = numFrames/(sim.tout(end)); % Adjust as needed
video.FrameRate = frameRate;

% Open the VideoWriter
open(video);

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];

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
    grid on
    title(['Flight time: ', num2str(sim.tout(i)), ' s'])

    xlim([min(R_inertial(:,1))-1 max(R_inertial(:,1))+1])
    ylim([min(R_inertial(:,2))-1 max(R_inertial(:,2))+1])
    zlim([min(R_inertial(:,3))-1 max(R_inertial(:,3))+1])
     % Capture the current frame
    frame = getframe(gcf);

    % Write the frame to the video
    writeVideo(video, frame);
end
% Close the VideoWriter
close(video);

%% Trajectory and attitude - Quaternions
% phi = sim.EulAng(1,:);
% theta = sim.EulAng(2,:);
% psi = sim.EulAng(3,:);
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
R_inertial = zeros(size(R_inertial));

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
for i = 1:10:length(sim.tout)
    % R1 = [
    % 1, 0, 0;
    % 0, cos(phi(i)), sin(phi(i));
    % 0, -sin(phi(i)), cos(phi(i))
    % ];
    % R2 = [
    % cos(theta(i)), 0, -sin(theta(i));
    % 0, 1, 0;
    % sin(theta(i)), 0, cos(theta(i))
    % ];
    % R3 = [
    % cos(psi(i)), sin(psi(i)), 0;
    % -sin(psi(i)), cos(psi(i)), 0;
    % 0, 0, 1
    % ];
    % T0 = R1*R2*R3;
    % invT0 = transpose(T0);
    invT0 = transpose(sim.DCM(:,:,i));
    attx = invTI0*invT0*x_v;
    atty = invTI0*invT0*y_v;
    attz = invTI0*invT0*z_v;

    figure(1);
    % plot3(R_inertial(:,1),R_inertial(:,2),R_inertial(:,3),'k-')
    % hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attx(1), attx(2), attx(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), atty(1), atty(2), atty(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attz(1), attz(2), attz(3), 'r')    
    xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
    hold off
    legend('Trajectory', 'x_b', 'y_b', 'z_b','interpreter', 'TeX')
    grid on
    % pause(.00000001)
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
%% Plot energies

figure(1)
plot(sim.P_E,'LineWidth',1.5)

hold on
plot(sim.K_E,'LineWidth', 1.5)

hold on
plot(sim.R_E,'LineWidth', 1.5)

title('Energy over time')
xlabel('Time [s]')
ylabel('Energy [J]')
legend('Potential Energy', 'Kinetic Energy', 'Rotational Energy')

%% IMU Measurements vs Ground Truth
close all

R_inertial = sim.R_inertial(:,:)';
U = sim.U(:,:)';
normU = vecnorm(U');
u = sim.u(:,:);
omega = sim.omega;
sample = 1/(f*dt);
t_s = sim.tout(1:sample:end);

acc = sim.acc(:,:);
R_s = R_inertial(1:sample:end,:);
u_s = u(:,1:sample:end);
acc_s = acc(:,1:sample:end);
omega_s = omega(1:sample:end,:);

acc_nav = sim.acc_nav;
omega_nav = sim.omega_nav;

u_nav = sim.u_nav;
u_nav = u_nav(1:sample:end, :);
R_nav = sim.R_nav(:,:);
R_nav = R_nav(:,1:sample:end);
EulAng_nav = sim.EulAng_nav(:,:);
EulAng = sim.EulAng;

err_vel = vecnorm(u_s - u_nav');
err_pos = vecnorm(R_s' - R_nav);
err_omega = vecnorm(omega_s' - omega_nav');

figure(1)
plot(sim.tout, acc(1,:), 'k-')
hold on
plot(t_s, acc_nav(:,1), 'k--')
legend('Ground truth', 'IMU estimates')
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
title('Acceleration along x - Actual vs estimated comparison')

figure(2)
plot(sim.tout, acc(2,:), 'k-')
hold on
plot(t_s, acc_nav(:,2), 'k--')
legend('Dynamic model data', 'IMU estimates')
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
title('Acceleration along y - Actual vs estimated comparison')

figure(3)
plot(sim.tout, acc(3,:), 'k-')
hold on
plot(t_s, acc_nav(:,3), 'k--')
legend('Dynamic model data', 'IMU estimates')
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
title('Acceleration along z - Actual vs estimated comparison')

figure(4)
plot(t_s, err_omega, 'k')
xlabel('Time[s]')
ylabel('Error [rad/s]')
title('Absolute Error in angular velocities estimation')

figure(5)
plot(sim.tout, normU, 'k-')
hold on
plot(t_s, vecnorm(u_nav'), 'k--')
legend('Ground truth', 'IMU estimates')
xlabel('Time [s]')
ylabel('Speed [m/s]')
title('Speed - Actual vs estimated comparison')

figure(6)
plot(R_inertial(:,1),R_inertial(:,2), 'k-')
hold on
plot(R_nav(1,:), R_nav(2,:), 'k--')
legend('Ground Truth', 'IMU estimates')
axis equal
xlabel('Displacemenent along x')
ylabel('Displacemenent along y')
title('Top view trajectory - Actual vs estimated comparison')

figure(7)
plot(t_s, err_vel)
xlabel('Time[s]')
ylabel('Error [m/s]')
title('Absolute Error in speed estimation')

figure(8)
plot(t_s, err_pos, 'k-')
xlabel('Time [s]')
ylabel('Error [m]')
title('Absolute Error in position estimation')