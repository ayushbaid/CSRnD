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

noise_sd = 1;



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

[f_omp, ~] = omp(y, A, noise_sd*sqrt(m), num_sinusoids*2);
[f_l1, ~] = l1solver(y, A, noise_sd*sqrt(m));
[f_iht, ~] = iht(y, A, noise_sd*sqrt(m), num_sinusoids*2);

x_omp = basis_mat*f_omp;
x_l1 = basis_mat*f_l1;
x_iht = basis_mat*f_iht;

%% Get L2 error in frequency domain
err_omp = norm(abs(f_original) - abs(f_omp))/norm(abs(f_original));
err_l1 = norm(abs(f_original) - abs(f_l1))/norm(abs(f_original));
err_iht = norm(abs(f_original) - abs(f_iht))/norm(abs(f_original));

fprintf('OMP relative error = %f\n', err_omp);
fprintf('L1 solver relative error = %f\n', err_l1);
fprintf('Iterative hard thresholding relative error = %f\n', err_iht);

%% Compare sparsity
sparse_omp = 1-sparsity_comp(f_omp, f_original, num_sinusoids*2)/(num_sinusoids*2);
sparse_l1 = 1-sparsity_comp(f_l1, f_original, num_sinusoids*2)/(num_sinusoids*2);
sparse_iht = 1-sparsity_comp(f_iht, f_original, num_sinusoids*2)/(num_sinusoids*2);

fprintf('OMP sparsity miss index = %f\n', sparse_omp);
fprintf('L1 solver sparsity miss index = %f\n', sparse_l1);
fprintf('IHT sparsity miss index = %f\n', sparse_iht);


%% Plots

figure();
subplot(231);
plot(abs(f_original), 'r');
legend('Original');
subplot(232);
plot(abs(f_noisy), 'k');
legend('Noisy');
subplot(233);
plot(abs(f_omp), 'b');
legend('OMP');
subplot(234);
plot(abs(f_l1), 'm');
legend('L1 solver');
subplot(235);
plot(abs(f_iht), 'c');
legend('IHT');


figure();
hold on;
plot(abs(f_original), 'r');
plot(abs(f_noisy), 'k');
plot(abs(f_omp), 'b');
plot(abs(f_l1), 'm');
plot(abs(f_iht), 'c');
hold off;
legend('Original', 'Noisy', 'OMP', 'L1 solver', 'IHT');
