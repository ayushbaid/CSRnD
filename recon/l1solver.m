function [f, r] = l1solver(y, A, e)
    % Solve the LP problem
    opts = spgSetParms('verbosity', 0);
    f = spg_bpdn(A, y, e, opts);

    r = y - A*f;
