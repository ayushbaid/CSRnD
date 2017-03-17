function [f, r] = l1_magic(y, A, e)
    addpath('../l1magic/Optimization');
    addpath('../l1magic/Data');

    % Do initial guess (min. energy)
    f = A\y;

    % Solve the LP problem
    f = l1qc_logbarrier(f, A, [], y, e);

    r = y - A*f;
