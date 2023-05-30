function [theta_hat_new] = evms_nSource_fn(X, theta_hat)
%% Separation and Tracking of Multiple Broadband Sources
% https://www.ese.wustl.edu/~nehorai/paper/ieeetae02.pdf

%% Generic sensor response:
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

%% Obtain Estimates of the incoming signal
% theta_hat = [pi/6 2*pi/3 5*pi/6; % azimuth                  [-pi, pi]
%         pi/6 2*pi/3 5*pi/6;      % elevation                [-pi/2,pi/2]
%         pi/6 2*pi/3 5*pi/6;      % polarization orientation [-pi/2,pi/2]
%         pi/6 2*pi/3 5*pi/6];     % polarization ellipticity [-pi/4,pi/4]
numSignals = 3;
S_theta_hat = zeros(6,3);
for i = 1:numSignals
    S_theta_hat(:,i) = double(subs(S_theta_n, [phi_n psi_n alpha_n beta_n], ...
         [theta_hat(1,i) theta_hat(2,i) theta_hat(3,i) theta_hat(4,i)]));
end

%% Form the pre-processing matrix from the estimates
% For now, assume that the desired source is theta_1 (repeat for 2 and 3)
[U, S, V] = svd(S_theta_hat);
% Use (47) to solve for k
syms k1;
P_H = zeros(3,4,6);
k = zeros(1,3);
for ind = 1:numSignals
    eq = abs(k1*V(ind,1)/S(1,1))^2 + abs(k1*V(ind,2)/S(2,2))^2 + ...
        abs(k1*V(ind,3)/S(3,3))^2 == 1;
    k(ind) = solve(eq,k1);
    % Combine the relavant equation to form the Preprocessing matrix
    P_H(ind,:,:) = [k(ind)*(1/S(1,1))*V(ind,1) k(ind)*(1/S(2,2))*V(ind,2) k(ind)*(1/S(3,3))*V(ind,3) 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1] * U'; % This is the hermitian transpose the the matrix
end
%% Calculate the Psuedo-inverses of the estimates
S_dot_theta = [
    gradient(S_theta_n(1),phi_n) gradient(S_theta_n(1),psi_n) gradient(S_theta_n(1),alpha_n) gradient(S_theta_n(1),beta_n)
    gradient(S_theta_n(2),phi_n) gradient(S_theta_n(2),psi_n) gradient(S_theta_n(2),alpha_n) gradient(S_theta_n(2),beta_n)
    gradient(S_theta_n(3),phi_n) gradient(S_theta_n(3),psi_n) gradient(S_theta_n(3),alpha_n) gradient(S_theta_n(3),beta_n)
    gradient(S_theta_n(4),phi_n) gradient(S_theta_n(4),psi_n) gradient(S_theta_n(4),alpha_n) gradient(S_theta_n(4),beta_n)
    gradient(S_theta_n(5),phi_n) gradient(S_theta_n(5),psi_n) gradient(S_theta_n(5),alpha_n) gradient(S_theta_n(5),beta_n)
    gradient(S_theta_n(6),phi_n) gradient(S_theta_n(6),psi_n) gradient(S_theta_n(6),alpha_n) gradient(S_theta_n(6),beta_n)
        ];
Q_plusR = zeros(3,4,6);
for ind = 1:numSignals
    Q = squeeze(P_H(ind,:,:)) * double(subs(S_dot_theta, [phi_n psi_n alpha_n beta_n], ...
             [theta_hat(1,ind) theta_hat(2,ind) theta_hat(3,ind) theta_hat(4,ind)]));
    Q_R = [real(Q(2,1)) real(Q(2,2)) real(Q(2,3)) real(Q(2,4));
               real(Q(3,1)) real(Q(3,2)) real(Q(3,3)) real(Q(3,4));
               real(Q(4,1)) real(Q(4,2)) real(Q(4,3)) real(Q(4,4));
               imag(Q(2,1)) imag(Q(2,2)) imag(Q(2,3)) imag(Q(2,4));
               imag(Q(3,1)) imag(Q(3,2)) imag(Q(3,3)) imag(Q(3,4));
               imag(Q(4,1)) imag(Q(4,2)) imag(Q(4,3)) imag(Q(4,4))];
    
    Q_plusR(ind,:,:) = inv(Q_R.' * Q_R) * Q_R.';
end
%% Find the optimal weight vector by solving the eigenvector problem
tLength = length(X(1,:));
R = X * X' / tLength;
W_P = zeros(4,3);
for i = 1:numSignals
    [N,V] = eig(squeeze(P_H(i,:,:)) * R * squeeze(P_H(i,:,:))','vector');
    % sort the eigenvalues
    [~,ind] = sort(V, 'descend');
    N = N(:,ind);
    % The principle (maximum) eigenvector is now the first eigenvector
    W_P(:,i) = N(:,1);
    
end
%% Get the new estimates
W_PR = zeros(6,3);
delta_theta = zeros(4,3);
theta_hat_new = zeros(4,3);
for ind = 1:numSignals
    W_PR(:,ind) = [real(W_P(2,ind) / W_P(1,ind));
             real(W_P(3,ind) / W_P(1,ind));
             real(W_P(4,ind) / W_P(1,ind));
             imag(W_P(2,ind) / W_P(1,ind));
             imag(W_P(3,ind) / W_P(1,ind));
             imag(W_P(4,ind) / W_P(1,ind))];
    
    delta_theta(:,ind) = squeeze(Q_plusR(ind,:,:)) * W_PR(:,ind);
    % Update the estimates
    theta_hat_new(:,ind) = theta_hat(:,ind) + delta_theta(:,ind);
    % roll into the valid range
    theta_hat_new(1,ind) = wrapToPi(theta_hat_new(1,ind));
    theta_hat_new(2,ind) = 2 * wrapToPi(theta_hat_new(2,ind)/2);
    theta_hat_new(3,ind) = 2 * wrapToPi(theta_hat_new(3,ind)/2);
    theta_hat_new(4,ind) = 4 * wrapToPi(theta_hat_new(4,ind)/4);
end
end