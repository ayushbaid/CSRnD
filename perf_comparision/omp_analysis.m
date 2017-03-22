%% Housekeeping

clear all;
close all;
clc;

addpath(genpath(pwd));

rng(0);


%% Define parameters
freq = [200, 400, 500, 600, 750];
coeffs = [1, 0.3, 2.3, 1.5, 0.8];
num_sinusoids = numel(freq);
n = 200;
fs = 2000;

solver = 'omp';

m_list = [1, 3, 10, 30, 50, 70, 90, 100];

noise_sd_list = [0, 3, 5, 10, 15, 30, 50, 100, 200];


%% Generate experiment params
num_exp = numel(m_list)*numel(noise_sd_list);

exp_params = zeros(num_exp, 2);
rrmse = zeros(num_exp, 1);

exp_idx = 1;
for idx1=1:numel(m_list)
    for idx2=1:numel(noise_sd_list);
        exp_params(exp_idx, 1) = m_list(idx1);
        exp_params(exp_idx, 2) = noise_sd_list(idx2);
        exp_idx = exp_idx + 1;
    end
end



%% Generate the ground truth
i = 0:(n-1);

x_original = zeros(n, 1);
for idx=1:num_sinusoids
    x_original = x_original + coeffs(idx)*sin(2*pi*freq(idx)*i/fs)';
end

sig_power= (sum(x_original.^2))/n;


%% Generate/Load the measurement matrix and the basis matrix

% Generate a Gaussian IID matrix
sensing_mat_full = randn(n, n);
basis_mat = dftmtx(n);
inv_basis_mat = conj(basis_mat)/n;


for exp_idx=1:num_exp
    fprintf('Exp #%d\n', exp_idx);

    m = ceil(exp_params(exp_idx, 1)*n/100);
    noise_sd = sqrt(exp_params(exp_idx, 2)*sig_power/100);

    x = x_original + randn(size(x_original))*noise_sd;

    sensing_mat = sensing_mat_full(1:m,:);
    A = sensing_mat*basis_mat;


    % Get the measurements
    y = sensing_mat*x;


    % Solve
    f_original = inv_basis_mat*x_original;
    f_noisy = inv_basis_mat*x;

    [f_omp, r] = omp(y, A, noise_sd*sqrt(m));

    x_omp = basis_mat*f_omp;

    % Get L2 error in frequency domain
    rrmse(exp_idx) = norm(abs(f_original) - abs(f_omp))/norm(abs(f_original));

    % fprintf('OMP relative error = %f\n', err_omp);

end


%% Plots
figure()
pcolor(m_list, noise_sd_list, reshape(rrmse, numel(noise_sd_list), []));
xlabel('Measured points (%)');
ylabel('Noise Power (%)');
title('OMP performance bechmark');
colorbar;

export_fig 'omp' '-dpng'


%% Save the rrmse
save('omp.mat', 'exp_params', 'rrmse');
