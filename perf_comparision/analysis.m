%% Housekeeping

clear all;
close all;
clc;

addpath(genpath(pwd));

rng('shuffle');


%% Define parameters
freq = [200, 400, 500, 600, 750];
coeffs = [1, 0.3, 2.3, 1.5, 0.8];
num_sinusoids = numel(freq);
n = 200;
fs = 2000;

% Top p1 peaks in original should be in top p2 of ouput
sparsity_comp_params = [2, 4];

m_list = [1, 3, 10, 30, 50, 70, 90, 100];

noise_sd_list = [0, 5, 15, 30, 60, 100, 150, 200];


%% Generate experiment params
num_reps = 10; % Number of times to repeat the whole experiment
num_exp = numel(m_list)*numel(noise_sd_list);

exp_params = zeros(num_exp, 2);
rrmse = zeros(num_exp, 3);
miss_index = zeros(num_exp, 3);

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

for rep_idx=1:num_reps
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

        [f_omp, ~] = omp(y, A, noise_sd*sqrt(m), num_sinusoids*2);
        [f_l1, ~] = l1solver(y, A, noise_sd*sqrt(m));
        [f_iht, ~] = iht(y, A, noise_sd*sqrt(m), num_sinusoids*2);

        % Get L2 error in frequency domain
        rrmse(exp_idx, 1) = rrmse(exp_idx, 1) + ...
                norm(abs(f_original) - abs(f_omp))/norm(abs(f_original));
        rrmse(exp_idx, 2) = rrmse(exp_idx, 2) + ...
                norm(abs(f_original) - abs(f_l1))/norm(abs(f_original));
        rrmse(exp_idx, 3) = rrmse(exp_idx, 3) + ...
                norm(abs(f_original) - abs(f_iht))/norm(abs(f_original));

        miss_index(exp_idx, 1) = miss_index(exp_idx, 1) + sparsity_comp(...
                f_omp, f_original, sparsity_comp_params(2), sparsity_comp_params(1));
        miss_index(exp_idx, 2) = miss_index(exp_idx, 2) + sparsity_comp(...
                f_l1, f_original, sparsity_comp_params(2), sparsity_comp_params(1));
        miss_index(exp_idx, 3) = miss_index(exp_idx, 3) + sparsity_comp(...
                f_iht, f_original, sparsity_comp_params(2), sparsity_comp_params(1));


    end
end

rrmse = rrmse./num_reps;
miss_index = miss_index./num_reps;


%% Plots
rrmse_max = max(rrmse(:));
figure('rend','painters','pos',[10 10 1500 1500])
subplot(221);
pcolor(m_list, noise_sd_list, reshape(rrmse(:, 1), numel(noise_sd_list), []));
xlabel('Measured points (%)');
ylabel('Noise Power (%)');
title('OMP');
caxis([0, rrmse_max]);
colorbar;
subplot(222);
pcolor(m_list, noise_sd_list, reshape(rrmse(:, 2), numel(noise_sd_list), []));
xlabel('Measured points (%)');
ylabel('Noise Power (%)');
title('L1 solver');
caxis([0, rrmse_max]);
colorbar;
subplot(223);
pcolor(m_list, noise_sd_list, reshape(rrmse(:, 3), numel(noise_sd_list), []));
xlabel('Measured points (%)');
ylabel('Noise Power (%)');
title('Iterative Hard Thres.');
caxis([0, rrmse_max]);
colorbar;


export_fig 'rrmse' '-dpng'

figure('rend','painters','pos',[10 10 1500 1500])
subplot(221);
pcolor(m_list, noise_sd_list, reshape(miss_index(:, 1), numel(noise_sd_list), []));
xlabel('Measured points (%)');
ylabel('Noise Power (%)');
title('OMP');
colorbar;
subplot(222);
pcolor(m_list, noise_sd_list, reshape(miss_index(:, 2), numel(noise_sd_list), []));
xlabel('Measured points (%)');
ylabel('Noise Power (%)');
title('L1 solver');
colorbar;
subplot(223);
pcolor(m_list, noise_sd_list, reshape(miss_index(:, 3), numel(noise_sd_list), []));
xlabel('Measured points (%)');
ylabel('Noise Power (%)');
title('Iterative Hard Thres.');
colorbar;
filename = sprintf('miss_%d_%d', sparsity_comp_params(1), ...
                   sparsity_comp_params(2));
export_fig(filename, '-dpng');

%% Save
save(sprintf('analysis_%d_%d.mat', sparsity_comp_params(1), ...
        sparsity_comp_params(2)), 'exp_params', 'rrmse', 'miss_index');
