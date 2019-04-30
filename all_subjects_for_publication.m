%% load data and calculate cohort connectome
% We created a matlab file for which all the connectomes have been read in,
% normalized (max weight=1 for each connectome) and then saved as a 3D array 
% (num_subjects, num_nodes, num_nodes)
% in addition there needs to be an age-vector for each subjects age
load('rockland.mat')

%rearange for easier indexing (num_nodes, num_nodes, num_subjects)
W = permute(W,[2,3,1]);

% define threshold in how many subjects edges need to be present
percentage = 0.9;

% find edges that exist in over 90% of the subjects
overall_group_connectome = sum(W(:,:,:)>0, 3)>(percentage*size(W(:,:,:),3));
% weight edges by average weight across cohort
overall_group_connectome = overall_group_connectome.*mean(W(:,:,:),3);

% load index that is false for brainstem and cerebellum and true otherwise
node_subset = csvread('remove_brains_cereb.csv')==1;

% cut down cohort connectome
overall_group_connectome = overall_group_connectome(node_subset,:);
overall_group_connectome = overall_group_connectome(:,node_subset);

%% calculate NDI
NDI_score = get_NDI(overall_group_connectome,max(overall_group_connectome(:)));

% define number of Gaussian for further analysis
num_gaussians = 3;

%% plot
figure

% create normalized histogram
[heights,locations] = hist(log(NDI_score(NDI_score>0)),40);
width = locations(2)-locations(1);
heights = heights / (sum(heights)*width);
h = bar(locations,heights,'hist');
set(h,'FaceColor',[0,0.5,0.75],'EdgeColor',[0, 0.5, 0.75])
xlabel('ln(NDI)')
ylabel('Probability')
hold on

% run gmm
gmfit = fitgmdist(log(NDI_score(NDI_score>0)),num_gaussians, 'Options', statset('MaxIter',1000));
mu = gmfit.mu;
sd = gmfit.Sigma;
pcomp = gmfit.PComponents;

% plot distribution
xx = linspace(-18, -6,100)';
plot(xx,pdf(gmdistribution(mu,sd,pcomp),xx),'-','Color',[0.8500, 0.3250, 0.0980],'LineWidth',3)
scatter(mu,pdf(gmdistribution(mu,sd,pcomp),mu),60,'dk','f')
yy=pdf(gmdistribution(mu,sd,pcomp),mu);
for ii = 1:length(mu)
    text(mu(ii)+0.1,yy(ii)+0.015,num2str(mu(ii),'%0.2f'));
end
xlim([-18,-3])

%% assign nodes to tiers
mu = sort(mu); % our study [-14.76;-10.75;-7.93]
NDI_labels = get_NDI_labels(NDI_score, mu(1:num_gaussians));

%% plot tier assignment, no left/right differentiation
colours = [[255 255 255];[86 180 233];[213 94 0]]/255.;
get_large_hitmap(NDI_labels,{'Tier 4', 'Tier 3', 'Tier 2', 'Tier 1'},'Labels_170_NKI_with_header_noLR.csv',colours,0)

%% topological assessment of NDI subnetworks
% calculate network measures for each subnetwork
network_measures = zeros(size(W,3),3,length(unique(NDI_labels)));
unique_labels = unique(NDI_labels);
for ii = 1:size(W,3)
    ii

    % extract subject connectome and remove brain stem and cerebellum
    connectome = W(:,:,ii);
    connectome = connectome(:,node_subset);
    connectome = connectome(node_subset,:);
    
    % calculate measures for each subnetwork defined by its label
    for label =1:length(unique(NDI_labels))
        idx = NDI_labels==unique_labels(label);

        % find subnetwork
        subnetwork = connectome(:,idx);
        subnetwork = subnetwork(idx,:);
        
        % calculate measures
        network_measures(ii,1,label) = transitivity_wu(subnetwork);
        network_measures(ii,2,label) = efficiency_wei(subnetwork);
        network_measures(ii,3,label) = assortativity_wei(subnetwork,0); 
    end    
end

%% plotting network measures
figure
colours = [[213 94 0];[230 159 0];[86 180 233];[0 114 178]]/255.;
m_names = {'T','E','a'};
plot_handles = [];
all_coefficients = [];
all_p_values = [];
for imeasure =1:size(network_measures,2)
    for label = 1:length(unique(NDI_labels))
        subplot(1,3,imeasure)
        hold on
        scatter(age, network_measures(:,imeasure,label),20,colours(label,:),'f', 'Marker','o');
        linm = fitlm(age,network_measures(:,imeasure,label),'linear');
        p = linm.Coefficients.Estimate([2,1]);
        all_coefficients = [all_coefficients p];
        all_p_values = [all_p_values linm.Coefficients.pValue];
        this_handle=plot(1:100,polyval(p,1:100),'Color',colours(label,:));
        plot_handles = [plot_handles this_handle];
        ylabel(m_names{imeasure})
        xlabel('Age')
    end
end

% create legend and put it to the top
[hL, hObj] = legend([plot_handles(1:max(NDI_labels(:)))],{'Tier 1','Tier 2','Tier 3','Tier 4'},'Orientation','horizontal');
newPosition = [0.5 0.95 0.05 0.05];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
% make legend lines thicker
hL=findobj(hObj,'type','line');
set(hL,'linewidth',2)   

%% rc analysis
% phi = rich_club_wu(overall_group_connectome);
% 
% % get random RC
% num_randoms = 1000;
% random = zeros(num_randoms,length(phi));
% for ii = 1:num_randoms
%     ii
%     Arand = null_model_und_sign(overall_group_connectome);
%     phi_rand = rich_club_wu(Arand);
%     random(ii,:) = phi_rand(1:length(phi));
% end
% 
% all_phi_norm_p = ones(size(random,2),1);
% for ii = 1:size(random,2)
%     idx = isfinite(random(:,ii));
%     [dummy,p] = ttest(phi(ii)./random(idx,ii),1,'Tail','right');
%     all_phi_norm_p(ii) = p;
% end
% 
% phi_rand = mean(random,1);
% 
% % normalized RC
% phi_norm = phi./phi_rand;
% 
% figure
% plot(phi_norm)
% hold on
% xx = 1:length(phi_norm);
% sig_level = 0.05/max(sum(connectome>0));
% scatter(xx(all_phi_norm_p(:)<sig_level), phi_norm(all_phi_norm_p(:)<sig_level),'f');
%
% % Result: k=36-47

%% find RC subnetwork assignments
idx = sum(overall_group_connectome~=0)>=47;
rc = find(idx);

feeder = [];
for ii = rc
    tmp_idx = find(overall_group_connectome(ii,:)>0);
    for jj = tmp_idx
        if ~any(rc==jj)
            feeder = [feeder jj] ;
        end
    end
end
feeder = unique(feeder);
sort_idx = rc;
sort_idx = [sort_idx feeder];
seeder = setdiff(1:size(overall_group_connectome,1),sort_idx);

rc_labels = 3*ones(170,1);
rc_labels(feeder) = 2;
rc_labels(seeder) = 1;

%% plot label assignment / comparison
all_labels = zeros(size(rc_labels,1),2*size(rc_labels,2));
all_labels(:,1) = rc_labels==3;

% find top NDI regions (same number as RC)
[aa,bb] = sort(NDI_score(:),'descend');
all_labels(bb(1:sum(rc_labels==3)),2) = 2;

% define group names
plot_group_names = {'RC','NDI'};
colours = [[255 255 255];[86 180 233];[213 94 0]]/255.;
get_hitmap(all_labels,plot_group_names,'Labels_170_NKI_with_header.csv',colours,1)

% plot full NDI and RC assignments 
get_large_hitmap(NDI_labels,{'Tier 4','Tier 3','Tier 2','Tier 1'},'Labels_170_NKI_with_header_noLR.csv',colours([1,3],:),0,0)
get_large_hitmap(rc_labels,{'S','F','RC'},'Labels_170_NKI_with_header_noLR.csv',colours([1,2],:),1,0)

%% topological assessment of RC subnetworks
% calculate network measures for each subnetwork
network_measures = zeros(size(W,3),3,length(unique(rc_labels)));
unique_labels = unique(rc_labels);
for ii = 1:size(W,3)
    ii

    % extract subject connectome and remove brain stem and cerebellum
    connectome = W(:,:,ii);
    connectome = connectome(:,node_subset);
    connectome = connectome(node_subset,:);
    
    % calculate measures for each subnetwork defined by its label
    for label =1:length(unique(rc_labels))
        idx = rc_labels==unique_labels(label);

        % find subnetwork
        subnetwork = connectome(:,idx);
        subnetwork = subnetwork(idx,:);
        
        % calculate measures
        network_measures(ii,1,label) = transitivity_wu(subnetwork);
        network_measures(ii,2,label) = efficiency_wei(subnetwork);
        network_measures(ii,3,label) = assortativity_wei(subnetwork,0); 
    end    
end

%% plotting network measures
figure
colours = [[0 114 178];[230 159 0];[213 94 0]]/255.;
m_names = {'T','E','a'};
plot_handles = [];
all_coefficients = [];
all_p_values = [];
for imeasure =1:size(network_measures,2)
    for label = 1:length(unique(rc_labels))
        subplot(1,3,imeasure)
        hold on
        scatter(age, network_measures(:,imeasure,label),20,colours(label,:),'f', 'Marker','o');
        linm = fitlm(age,network_measures(:,imeasure,label),'linear');
        p = linm.Coefficients.Estimate([2,1]);
        all_coefficients = [all_coefficients p];
        all_p_values = [all_p_values linm.Coefficients.pValue];
        this_handle=plot(1:100,polyval(p,1:100),'Color',colours(label,:));
        plot_handles = [plot_handles this_handle];
        ylabel(m_names{imeasure})
        xlabel('Age')
    end
end

% create legend and put it to the top
[hL, hObj] = legend([plot_handles(3:-1:1)],{'RC','F','S'},'Orientation','horizontal');
newPosition = [0.5 0.95 0.05 0.05];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
% make legend lines thicker
hL=findobj(hObj,'type','line');
set(hL,'linewidth',2)   

%% Analyse NDI tier assignment on subject level
% cut down matrix to exclude brain stem and cerebellum
sub_W = W(node_subset,:,:);
sub_W = sub_W(:,node_subset,:);

% create output variable
label_sub = zeros(sum(node_subset),size(sub_W,3));

% loop through subjects
for isub = 1:size(label_sub,2)
    isub
    % calculate NDI score for subject
    tmp_W = sub_W(:,:,isub);
    NDI_score = get_NDI(tmp_W,max(tmp_W(:)));

    % run gmm
    try
        gmfit = fitgmdist(log(NDI_score(NDI_score>0)),num_gaussians, 'Options', statset('MaxIter',1000));
    catch
        %do nothing
    end
    mu = gmfit.mu;
    sd = gmfit.Sigma;
    pcomp = gmfit.PComponents;
    mu = sort(mu);
    
    % calculate labels
    label_sub(:,isub) = get_NDI_labels(NDI_score, mu(1:num_gaussians));
end

%% plot subject level analysis

figure
[dummy, idx] = sort(median(label_sub, 2));
sorted_labels = label_sub(idx,:);

% add image of all nodal assignments
subplot(2,1,1)
[dummy, age_idx] = sort(age);
sorted_labels = sorted_labels(:,age_idx)';
imagesc(sorted_labels)
ylabel({'Subject index' '(from young to old)'})
xlabel('')
set(gca,'XTickLabel','')
cbh = colormap([colours; [0 0 0]]);
caxis([0.5,4.5])
h=colorbar;
ylabel(h,'Tier')
set(gca,'Position',[0.13, 0.53, 0.775, 0.4])

% add median plot
subplot(2,1,2)
errorbar(median(label_sub(idx,:), 2), std(label_sub(idx,:),[],2), 'LineStyle','none')
xlabel('# regions sorted by median assignment')
ylabel({'Median NDI' 'assignment'})
yticks([0:5])
grid on 
set(gca, 'YGrid', 'on', 'XGrid', 'off')
xlim([1 170])
set(gca,'Position',[0.13, 0.12, 0.775, 0.4])
set(gca, 'YTickLabel',{'','1','2','3','4',''})