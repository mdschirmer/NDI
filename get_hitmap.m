function get_hitmap(region_idx,group_names,file,colours,xlocation)

    % predefine colormap
    if colours==0
       colours = parula(length(unique(region_idx(:)))); 
    end

    labels = [];
    for group = 1:size(region_idx,2)
        labels = [labels; find(region_idx(:,group))];
    end

    % read in region namesregion_names
    data = readtable(file);
   
    %find unique labels
    unique_labels = sort(unique(labels));
    region_names = table2cell(data(unique_labels,2));
    
    % cutdown matrix
    region_idx_cutdown = region_idx(unique_labels,:);

    % sort regions (due to multiple nodes)
    [region_names, idx] = sort(region_names);
    region_idx_cutdown = region_idx_cutdown(idx, :);
    
    % sort by ignoring left/right label (alphabetic order)
    tmp={};
    for ii = 1:length(region_names)
        jj = strsplit(region_names{ii}, ' ');
        tmp(ii) = jj(2);
    end
   
    [dummy, idx] = sort(tmp);
    region_names = region_names(idx);
    region_idx_cutdown = region_idx_cutdown(idx,:);
    
    % find duplicates
    idx = zeros(size(region_names));
    for ii = 1:length(region_names)
        idx(ii) = find(strcmp(region_names,region_names(ii)),1);
    end 

    % create image to display
    unique_regions = unique(idx);
    hitmap = zeros(length(unique_regions),size(region_idx,2));
    region_names = region_names(unique_regions);
    
    for ii =1:length(unique_regions)
        hitmap(ii,:) = max(region_idx_cutdown(idx==unique_regions(ii),:),[],1);
    end
    
    %plot
    figure
    imagesc(hitmap')
    set(gca, 'YTick', 1:size(region_idx,2) ,'YTickLabel',group_names,'XTick',1:length(region_names),'XTickLabel',region_names,'XTickLabelRotation',-45)
    colormap(colours)
    
    % set x axis to be on top if requested
    if xlocation==1
        set(gca,'XAxisLocation','top')
    end
    
end