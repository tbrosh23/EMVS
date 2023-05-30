%% Separation and Tracking of Multiple Broadband Sources
% https://www.ese.wustl.edu/~nehorai/paper/ieeetae02.pdf
close all;
clear all;
sympref('FloatingPointOutput',true);
rng('default');

%% Sensor response:
syms phi_n psi_n alpha_n beta_n;
%  First 3 components are E field resonse, last 3 are B field response
S_theta_n = [-sin(phi_n)            -cos(phi_n)*sin(psi_n);
             cos(phi_n)             -sin(phi_n)*sin(psi_n);
             0                      cos(psi_n);
             -cos(phi_n)*sin(psi_n) sin(phi_n);
             -sin(phi_n)*sin(psi_n) -cos(phi_n);
             cos(psi_n)             0] * ...
            [cos(alpha_n)  sin(alpha_n);
             -sin(alpha_n) cos(alpha_n)] * ...
            [cos(beta_n);
             1i*sin(beta_n)];

%% Broadband signals received
numSignals = 3;
t = 0:.01:6.99;
a = zeros(numSignals, length(t));
omega1 = 70;
omega2 = 110;
omega3 = 50;
a(1,:) = 10^(23/20)*exp(1i * omega1 * (t-3*pi/4)); % convert from dB
a(2,:) = 10^(20/20)*exp(1i * omega2 * (t-pi/8));
a(3,:) = 10^(21/20)*exp(1i * omega3 * (t+2*pi/3));

theta_in = [30 100 -140; % azimuth
        15 80 -55; % elevation
        5  35 65; % polarization orientation
        -2 15 -12] * pi/180; % polarization ellipticity

plot3(t,real(a(1,:)), imag(a(1,:)))
hold on
plot3(t,real(a(2,:)), imag(a(2,:)))
plot3(t,real(a(3,:)), imag(a(3,:)))
hold off
S_theta = zeros(6,3);
for i = 1:numSignals
    S_theta(:,i) = double(subs(S_theta_n, [phi_n psi_n alpha_n beta_n], ...
             [theta_in(1,i) theta_in(2,i) theta_in(3,i) theta_in(4,i)]));
end
n = 10^(1/20)*(rand(6,length(t))-.5);
n = 10^(1/20)*wgn(6, length(t), 1);
X_t = a(1,:).*S_theta(:,1) + a(2,:).*S_theta(:,2) + a(3,:).*S_theta(:,3) + n;
figure
tiledlayout(3,2)
for i = 1:6
    nexttile
    plot3(t,real(X_t(i,:)), imag(X_t(i,:)))
end

%% Obtain Estimates of the incoming signal
theta_hat = [40 110 -140; % azimuth                  [-pi, pi]
        20 70 -55;      % elevation                [-pi/2,pi/2]
        -10 25 65;      % polarization orientation [-pi/2,pi/2]
        -5 20 -12]...   % polarization ellipticity [-pi/4,pi/4]
        * pi/180;    

%% Find new estimates
covar_length = 100;
theta_hat_new = zeros(4,3,floor(length(X_t(1,:))/covar_length)+1);
theta_hat_new(:,:,1) = theta_hat;
%theta_hat_new(:,:,1) = ones(4,3);
for window = 1:floor(length(X_t(1,:))/covar_length)
    theta_hat_new(:,:,window+1) = ...
        evms_nSource_fn(X_t(:,(covar_length*(window-1)+1):(covar_length*window)),...
            theta_hat_new(:,:,window));
    disp(['estimate ',window,': ']);
    theta_hat_new(:,:,window+1) * 180/pi
end

figure
tiledlayout(3,1);
for i = 1:3
    nexttile
    plot(1:floor(length(X_t(1,:))/covar_length)+1,squeeze(theta_hat_new(:,i,:)).');
    hold on
    ylim([-pi pi])
    legend(["azimuth","elevation","orientation","ellipticity"])
    hold off
end

