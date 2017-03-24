function count = sparsity_comp(f, f_ref, s)
%sparsity_comp Compares the peak locations of f to f_ref; counts the peaks of
% f_ref captured by f
% Input Arguments:
% * f: the query frequency response
% * f_ref: reference frequency response
% * s: the number of peaks to look for

[~, q_idx] = sort(abs(f), 'descend');
[~, ref_idx] = sort(abs(f_ref), 'descend');

q_peaks = q_idx(1:s);
ref_peaks = ref_idx(1:s);


% Count the common indices
indicator = false(length(f), 1);
indicator(ref_peaks) = true;

count = 0;

for i=1:s
    if indicator(q_peaks(i))
        count = count+1;
    end
end
