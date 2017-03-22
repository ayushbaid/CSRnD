%% Housekeeping

clear all;
close all;
clc;

rng(0);

addpath(genpath(pwd));


%% Define parameters
freq = [200, 400, 500, 600, 750];
coeffs = [1, 0.3, 2.3, 1.5, 0.8];
num_sinusoids = numel(freq);

n = 200;
m = 50;

noise_sd = 0.15;



%% Generate the signal
fs = 2000;
i = 0:(n-1);

x_original = zeros(n, 1);
for idx=1:num_sinusoids
    x_original = x_original + coeffs(idx)*sin(2*pi*freq(idx)*i/fs)';
end
x = x_original + randn(size(x_original))*noise_sd;

%% Generate/Load the measurement matrix and the basis matrix

% Generate a Gaussian IID matrix
sensing_mat = randn(m, n);
basis_mat = dftmtx(n);
A = sensing_mat*basis_mat;


%% Get the measurements
y = sensing_mat*x;


%% Solve
inv_basis_mat = conj(basis_mat)/n;
f_original = inv_basis_mat*x_original;
f_noisy = inv_basis_mat*x;

[f_omp, ~] = omp(y, A, noise_sd*sqrt(m));
[f_l1, ~] = l1solver(y, A, noise_sd*sqrt(m));

x_omp = basis_mat*f_omp;
x_l1 = basis_mat*f_l1;

%% Get L2 error in frequency domain
err_omp = norm(abs(f_original) - abs(f_omp))/norm(abs(f_original));
err_l1 = norm(abs(f_original) - abs(f_l1))/norm(abs(f_original));

fprintf('OMP relative error = %f\n', err_omp);
fprintf('L1 solver relative error = %f\n', err_l1);

%% Plots

figure();
subplot(221);
plot(abs(f_original), 'r');
legend('Original');
subplot(222);
plot(abs(f_noisy), 'k');
legend('Noisy');
subplot(223);
plot(abs(f_omp), 'b');
legend('OMP');
subplot(224);
plot(abs(f_l1), 'm');
legend('L1-magic');


figure();
hold on;
plot(abs(f_original), 'r');
plot(abs(f_noisy), 'k');
plot(abs(f_omp), 'b');
plot(abs(f_l1), 'm');
hold off;
legend('Original', 'Noisy', 'OMP', 'L1-magic');
