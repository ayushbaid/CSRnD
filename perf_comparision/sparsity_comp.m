function miss_index = sparsity_comp(f, f_ref, s_bucket, s_ref)
%sparsity_comp Compares the peak locations of f to f_ref; counts the peaks of
% f_ref captured by f
% Input Arguments:
% * f: the query frequency response
% * f_ref: reference frequency response
% * s_bucket: the number of peaks to search for in f
% * s_ref: the number of peaks of f_ref we want to recover

% As we have only real data, we will trim frequenciess
half_idx = length(f)/2;

[~, q_idx] = sort(abs(f(1:half_idx)), 'descend');
[~, ref_idx] = sort(abs(f_ref(1:half_idx)), 'descend');

q_peaks = q_idx(1:s_bucket);
ref_peaks = ref_idx(1:s_ref);


% Count the common indices
indicator = false(half_idx, 1);
indicator(q_peaks) = true;

count = 0;

for i=1:s_ref
    if indicator(ref_peaks(i))
        count = count+1;
    end
end

miss_index = 1-count/s_ref;
