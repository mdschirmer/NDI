function NDI = get_NDI(W, max_I)

    % calculate information measure for graph
    DG = distance_wei(1./W);
    % define maximum information measure based on minimum path length
    IG = 1./DG;
    % normalize
    IG = IG/max_I;

    % calculate path dependence for each knock-out node n
    DI_i = zeros(size(W));
    nodes = 1:size(W,1);
    for n = nodes
        % knock out node n
        reduced_idx = nodes~=n;
        Wreduced = W(reduced_idx,:);
        Wreduced = Wreduced(:, reduced_idx);

        % calcualte information measure for graph-n
        D = distance_wei(1./Wreduced);
        A = isfinite(D);
        I = 1./D;
        % normalize
        I = I/max_I;

        IGreduced = IG(reduced_idx,:);
        IGreduced = IGreduced(:,reduced_idx);

        % node pair measure
        DIij = IGreduced-I; 
        DIij(A~=1) = 1;
        DIij(1:1+size(DIij,1):end) = 0;

        % nodal measure
        nodedependencyidx = sum(DIij,2) / (size(W,1)-1.);
        DI_i(reduced_idx,n) = nodedependencyidx;

    end
    % nodal values as average of NDI
    NDI = sum(DI_i, 1) / (size(W,1)-1.);
    NDI = NDI';
    
end