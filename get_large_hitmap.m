function get_large_hitmap(region_idx,group_names,file,colours,xlocation, alphabetical)

    % check if alphabetical ordering was set
    if nargin<6
        alphabetical=false;
    end

    % need to flip labels in case the region names are supposed to be on
    % the other side
    if xlocation==1
        region_idx = max(region_idx) + 1 - region_idx;
    end
    
    % predefine colormap
    if colours==0
       colours = parula(length(unique(region_idx(:)))); 
    end

    % find the region indices for each subnetwork
    labels = [];
    for group = 1:size(region_idx,2)
        labels = [labels; find(region_idx(:,group))];
    end

    % read in region names
    data = readtable(file);
   
    %find unique labels
    unique_labels = sort(unique(labels));
    region_names = table2cell(data(unique_labels,2));

    % unique region names
    unique_regions = sort(unique(region_names));
    unique_labels = zeros(length(unique_regions), length(unique(region_idx)));
    for ii = 1:length(region_names)
        region_found = -1;
        jj=1;
        while region_found==-1
            if strcmp(unique_regions{jj},region_names{ii})
                region_found = 1;
            else
                jj=jj+1;
            end
        end
        unique_labels(jj,region_idx(ii)) = 1;
    end
    
    region_idx = unique_labels;
    
    [~, aaa] = sort(sum(unique_labels,2));
    unique_labels = unique_labels(aaa,:);
    unique_regions = unique_regions(aaa);

    % find number of regions that only appear in one group
    subset = 1:sum(sum(unique_labels,2)==1);
    
    swap_idx = [];
    for group = 1:size(unique_labels,2)
        swap_idx = [swap_idx; find(unique_labels(subset,group))];
    end
    swap_idx = [swap_idx; ((length(swap_idx) + 1):size(unique_labels,1))'];
    
    unique_labels = unique_labels(swap_idx,:);
    unique_regions = unique_regions(swap_idx);
    
    if alphabetical
       [unique_regions, idx] = sort(unique_regions);
       unique_labels = unique_labels(idx,:);
    end
    
    %plot
    figure 
    % set x axis to be on top if requested
    if xlocation==1
        imagesc(flipud(unique_labels'))
        set(gca, 'YTick', 1:size(region_idx,2) ,'YTickLabel',group_names,'XTick',1:length(unique_regions),'XTickLabel',unique_regions,'XTickLabelRotation',90, 'YTickLabelRotation',135)
        colormap(colours)
        set(gca,'XAxisLocation','top')
    else
        imagesc(flipud(unique_labels'))
        set(gca, 'YTick', (1:size(region_idx,2)) ,'YTickLabel',group_names,'XTick',1:length(unique_regions),'XTickLabel',unique_regions,'XTickLabelRotation',90, 'YTickLabelRotation',135)
        colormap(colours)
    end
    
end