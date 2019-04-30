function labels = get_NDI_labels(NDI_score, mu)

% make sure mu increases
mu = sort(mu);

% predetermined cutoffs
% cutoffs = [-inf (-14.76-10.75)/2. (-10.75-7.93)/2.  0];

% cutoffs
diff_mu = diff(mu);
cutoffs = [-inf; mu(1:length(diff_mu))+diff_mu; 0];

% initialize output
labels = zeros(size(NDI_score));

% set NDI=0 to label 1
labels(NDI_score==0) = 1;

% assign remaining labels based on cutoffs
for jj=2:(length(mu)+1)
    idx = log(NDI_score)>cutoffs(jj-1) & log(NDI_score)<=cutoffs(jj);
    labels(idx) = jj;
end

% change order so most important label is tier 1
labels = max(labels) + 1 - labels;

end