function [f, r] = iht(y, A, e, s)
    f = zeros(size(A,2), 1);
    r = y - A*f;
    norm_r = norm(r, 2);

    A_t = A';

    step_size = 1;
    step_size_threshold = 1e-10;

    while (norm_r > e && step_size>step_size_threshold)
        temp = f + step_size*A_t*r;

        % Perform hard thresholding
        [~, sort_idx] = sort(abs(temp), 'descend');
        select_idx = sort_idx(1:s);
        f_new = zeros(size(f));
        f_new(select_idx) = temp(select_idx);

        r_new = y - A*f_new;
        norm_r_new = norm(r_new, 2);

        if (norm_r_new < norm_r)
            f = f_new;
            r = r_new;
            norm_r = norm_r_new;

            step_size = step_size*1.2;
        else
            step_size = step_size*0.5;
        end
    end
