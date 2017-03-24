function [f,r] = omp(y, A, e, s)

m = size(A,1);
n = size(A,2);

% Initializations
r = y;
support_set = false(n, 1);

iter_idx = 0;

A_norms = sqrt(sum(A.^2, 1));

f = zeros(n,1);

while (iter_idx < s && norm(r, 2)>e)

%     projections = abs(bsxfun(@ldivide, r.'*A, A_norms));
    projections = abs(r'*A);

    [~, idx] = max(projections);
    support_set(idx) = true;

    A_ = zeros(size(A));
    A_(:, support_set) = A(:, support_set);

    f = pinv(A_)*y;
    r_new = y - A_*f;

    iter_idx = iter_idx + 1;

%     if (r==r_new)
%         break;
%     end

    r = r_new;
end
