%% load data and calculate cohort connectome
% We created a matlab file for which all the connectomes have been read in,
% normalized (max weight=1 for each connectome) and then saved as a 3D array 
% (num_subjects, num_nodes, num_nodes)
% in addition there needs to be an age-vector for each subjects age
load('rockland.mat')

%rearange for easier indexing (num_nodes, num_nodes, num_subjects)
W = permute(W,[2,3,1]);

% load index to remove brainstem and cerebellum 
node_subset = csvread('remove_brains_cereb.csv')==1;
W = W(node_subset,:,:);
W = W(:,node_subset,:);

%% create group connectomes
group_connectome = zeros(size(W,1),size(W,2),4);
% define age bounds for groups
groups = [0,20,40,60,inf];
% define threshold in how many subjects edges need to be present
percentage = 0.9;
% create output for subject's group label
agegroup_labels = zeros(size(W,3),1);

for ii = 1:(length(groups)-1)
    % find subjects with age within group
    idx = age>groups(ii) & age<=groups(ii+1);
    agegroup_labels(idx) = ii;
    
    % get connectomes within group
    subset = W(:,:,idx);
    nsubs = sum(idx);
    
    % find edges that exist in over 90% of the subjects within group
    Ag = sum(subset>0,3)>(percentage*nsubs);
    % weight edges by average weight across cohort
    group_connectome(:,:,ii) = Ag.*mean(subset,3);
end

% define group names
group_names = {'U20','U40','U60','O60'};

%% find number of Gaussians for GMM fit
% initialize output
NDI_scores = zeros(size(group_connectome,1),size(group_connectome,3));

% range of number of Gaussians to look at
g = 1:4;
nK = numel(g);

% initialize output
gm = cell(nK,size(group_connectome,3));
aic = zeros(nK,size(group_connectome,3));
bic = zeros(nK,size(group_connectome,3));

for ig = 1:size(group_connectome,3)
    ig
    % get group connectome
    groupW =  group_connectome(:,:,ig);
    
    % calculate NDI and store for later
    this_NDI_score = get_NDI(groupW,max(groupW(:)));
    NDI_scores(:,ig) = this_NDI_score;

    % estimate Gaussians
    for ii = 1:nK
        ii
        % estimate GMM model
        try
            gm = fitgmdist(log(this_NDI_score(this_NDI_score>0)),g(ii),'Options',statset('MaxIter',1000));
        catch
            % do nothing
        end
            
        % extract Akaike and Bayesian Information Criteria
        aic(ii, ig) = gm.AIC;
        bic(ii, ig) = gm.BIC;
    end
end

%% plot GMM g assessment
figure
colours = [[255 255 255];[86 180 233];[213 94 0]]/255.;

% create bar plot
hb = bar([mean(aic,1)' mean(bic,1)']);
set(hb(1), 'FaceColor',colours(2,:))
set(hb(2), 'FaceColor',colours(3,:))
hold on

% create error bars
xData = hb(1).XData+hb(1).XOffset;
errorbar(xData,mean(aic,1),std(aic,1),'LineStyle','none','Color','k')
xData = hb(2).XData+hb(2).XOffset;
errorbar(xData,mean(bic,1),std(bic,1),'LineStyle','none','Color','k')

title('Evaluation for various number of mixtures $g$','Interpreter','latex');
xlabel('$g$','Interpreter','Latex');
ylabel('Information');
set(gca,'xticklabel',g)
grid on
ylim([min([bic(:);aic(:)])-1, max([bic(:);aic(:)])+1])
legend('AIC', 'BIC')

%% determine number of Gaussians
num_gaussians = 3;

%% Calculate NDI tier assginments and plot distributions
% initialize output
NDI_labels = zeros(size(group_connectome,1),size(group_connectome,3));
label_mu = zeros(num_gaussians,size(group_connectome,3));
% our study: label_mu = [[-8.29; -11.25;-14.34] [-8.06; -11.10;-15.82] [-8.06; -11.06; -14.73] [-8.33; -11.16;; -14.51]]

figure
for group =1:size(group_connectome,3)
    group
    NDI_score = NDI_scores(:,group);
    
    % GMM fit
    gmfit = fitgmdist(log(NDI_score(NDI_score>0)),num_gaussians, 'Options', statset('MaxIter',1000));
    label_mu(:,group) = gmfit.mu(:);
    sd = gmfit.Sigma;
    pcomp = gmfit.PComponents;

    % plot
    subplot(2,2,group)
    hold on
    % normalized histogram
    [heights,locations] = hist(log(NDI_score(NDI_score>0)),30);
    width = locations(2)-locations(1);
    heights = heights / (sum(heights)*width);
    h = bar(locations,heights,'hist');
    set(h,'FaceColor',[0,0.5,0.75],'EdgeColor',[0, 0.5, 0.75])
    xlabel('ln(NDI)')
    ylabel('Probability')
    title(group_names{group})
    
    % GMM distributions
    xx = linspace(min(locations), max(locations),100)';
    plot(xx,pdf(gmdistribution(label_mu(:,group),sd,pcomp),xx),'-','Color',[0.8500, 0.3250, 0.0980],'LineWidth',3)
    scatter(label_mu(:,group),pdf(gmdistribution(label_mu(:,group),sd,pcomp),label_mu(:,group)),60,'dk','f')
    
    % add labels to peaks
    yy=pdf(gmdistribution(label_mu(:,group),sd,pcomp),label_mu(:,group));
    for ii = 1:length(label_mu(:,group))
        text(label_mu(ii,group)-2,yy(ii)+0.05,num2str(label_mu(ii,group),'%0.2f'));
    end
    
    % set limits
    xlim([-18,-6])
    ylim([0, 0.4])

    % assign labels
    NDI_labels(:,group) = get_NDI_labels(NDI_score, label_mu(:,group));
end

%% topological assessment of NDI subnetworks
% set output
network_measures = zeros(size(W,3),3,num_gaussians+1);

% loop through groups
for group =1:length(group_names)
    group
    labels = NDI_labels(:,group);
    unique_labels = unique(labels);
    
    % find subjects within group
    for ii = find(age>groups(group) & age<=groups(group+1))
        ii
        % get subject's connectome
        connectome = W(:,:,ii);
        
        % calculate for each subnetwork defined by its label
        for label =1:length(unique(labels))
            idx = labels==unique_labels(label);
            
            % find subnetwork
            subnetwork = connectome(:,idx);
            subnetwork = subnetwork(idx,:);
            
            % calculate measures
            network_measures(ii,1,label) = transitivity_wu(subnetwork);
            network_measures(ii,2,label) = efficiency_wei(subnetwork);
            network_measures(ii,3,label) = assortativity_wei(subnetwork,0);
        end
    end
end

%% plotting
figure
hold all
m_names = {'T','E','a'};
labels = unique(NDI_labels(:));

% set titles
for ii = 1:4
    subplot(3,length(unique(labels)),ii)
    hold on
    title(['Tier ' num2str(ii)])
end

% set measure names
for ii = 1:3
    subplot(3,length(unique(labels)),(ii-1)*length(unique(labels))+1)
    hold on
    ylabel(m_names{ii})
end

% adjust subplot size
subplot_width=0.18;
subplot_height=0.25;
subplot_x = [0.1 0.3 0.5 0.7];
subplot_y = [0.70 0.43 0.16];

% plot each measure
for imeasure =1:3
    % for each tier
    for label = 1:length(unique(labels))
        subplot(3,length(unique(labels)),label+(imeasure-1)*length(unique(labels)))
        hold on
        boxplot(network_measures(:,imeasure,label),agegroup_labels, 'Labels',group_names);
        ylim([min(min(network_measures(:,imeasure,:))) max(max(network_measures(:,imeasure,:)))])
        
        % remove unwanted x and y labels
        if (label+(imeasure-1)*4) ~= 1 && (label+(imeasure-1)*4) ~= 5 && (label+(imeasure-1)*4) ~= 9
            set(gca, 'YTickLabel','')
        end
        if (label+(imeasure-1)*4) < 9
            set(gca, 'XTickLabel','')
        end
        
        % adjust subplot size
        set(gca,'Position',[subplot_x(label), subplot_y(imeasure), subplot_width, subplot_height])
    end
end


%% Run RC analysis
% figure
% for group =1:size(group_connectome,3)
%     subplot(2,2,group)
%     group
%     
%     % get group connectome
%     connectome = group_connectome(:,:,group);
%     
%     % calculate RC
%     phi = rich_club_wu(connectome);
%     
%     % get random RC
%     nrand=1000;
%     random = inf*ones(nrand,length(phi));
%     for ii = 1:nrand
%         ii
%         Arand = null_model_und_sign(connectome);
%         phi_rand = rich_club_wu(Arand);
%         random(ii,:) = phi_rand;
%     end
%     
%     % calculate significance
%     all_phi_norm_p = ones(size(random,2),size(group_connectome,3));
%     for ii = 1:size(random,2)
%         idx = isfinite(random(:,ii));
%         if sum(idx) ~= 0 || isfinite(phi(ii))
%             [h,p] = ttest(phi(ii)./random(idx,ii),1,'Tail','right');
%             all_phi_norm_p(ii,group) = p;
%         end
%     end
%     
%     % get average random phi
%     phi_rand = mean(random,1);
%     
%     % calculate normalized RC and plot
%     phi_norm = phi./phi_rand;
%     plot(phi_norm)
%     hold on
%     
%     % scatter dots for significant values
%     xx = 1:length(phi_norm);
%     sig_level = 0.05/max(sum(connectome>0));
%     scatter(xx(all_phi_norm_p(:,group)<sig_level), phi_norm(all_phi_norm_p(:,group)<sig_level),'f');
%     
% end

% % Result:
% % U20: k=34-49
% % U40: k=36-48
% % U60: k=37-47
% % O60: k=39-53
% %
% 
% kU20 = 34:49;
% kU40 = 36:48;
% kU60 = 37:47;
% kO60 = 39:53;
% rc_degrees = {kU20, kU40, kU60, kO60};
% 
% figure
% hold on
% for ii = 1:length(rc_degrees)
%    
%     x = [];
%     y = [];
%     
%     kThisGroup = rc_degrees{ii};
%     
%     for jj = 1:length(kThisGroup)
%        y = [y; kThisGroup(jj)-(-1)^ii * 0.1];
%        x = [x; sum(sum(group_connectome(:,:,ii)~=0)>=kThisGroup(jj))-(-1)^ii * 0.1];
%     end
%     
%     scatter(x,y,'f')
%     xlabel('Size of RC')
%     ylabel('k for Group')
%     
% end
% plot([15, 15], [34, 53],'k')
% legend({'U20','U40', 'U60', 'O60'})
% ylim([33.9, 53.1])
% xlim([4.9 25.1])
% grid minor
%
% % all groups have RC of size 15 at  
% kmax =[44 45 46 44];

kmax =[49 48 47 53];

%% find RC subnetwork assignments
rc_labels = zeros(size(W,1),length(kmax));
for group =1:length(kmax)
    
    % find RC
    idx = sum(group_connectome(:,:,group)~=0)>=kmax(group);
    rc = find(idx);

    feeder = [];
    for ii = rc
        tmp_idx = find(group_connectome(ii,:,group)>0);
        for jj = tmp_idx
            if ~any(rc==jj)
                feeder = [feeder jj] ;
            end
        end
    end
    feeder = unique(feeder);
    sort_idx = rc;
    sort_idx = [sort_idx feeder];
    seeder = setdiff(1:size(group_connectome,1),sort_idx);

    rc_labels(rc, group) = 3;
    rc_labels(feeder, group) = 2;
    rc_labels(seeder,group) = 1;
end

%% topological assessment of NDI subnetworks
% set output
network_measures = zeros(size(W,3),3,num_gaussians+1);

% loop through groups
for group =1:length(group_names)
    group
    labels = rc_labels(:,group);
    unique_labels = unique(labels);
    
    % find subjects within group
    for ii = find(age>groups(group) & age<=groups(group+1))
        ii
        % get subject's connectome
        connectome = W(:,:,ii);
        
        % calculate for each subnetwork defined by its label
        for label =1:length(unique(labels))
            idx = labels==unique_labels(label);
            
            % find subnetwork
            subnetwork = connectome(:,idx);
            subnetwork = subnetwork(idx,:);
            
            % calculate measures
            network_measures(ii,1,label) = transitivity_wu(subnetwork);
            network_measures(ii,2,label) = efficiency_wei(subnetwork);
            network_measures(ii,3,label) = assortativity_wei(subnetwork,0);
        end
    end
end

%% plotting
figure
hold all
m_names = {'T','E','a'};
subnetwork_names = {'RC','F','S'};
% set titles
for ii = 1:3
    subplot(3,length(unique(rc_labels)),ii)
    hold on
    title(subnetwork_names{ii})
end

% set measure names
for ii = 1:3
    subplot(3,length(unique(rc_labels)),(ii-1)*length(unique(rc_labels))+1)
    hold on
    ylabel(m_names{ii})
end

% adjust subplot size
subplot_width=0.18;
subplot_height=0.25;
subplot_x = [0.1 0.3 0.5 0.7];
subplot_y = [0.70 0.43 0.16];

% plot each measure
for imeasure =1:3
    % for each tier
    for label = 1:length(unique(rc_labels))
        subplot(3,length(unique(rc_labels)),label+(imeasure-1)*length(unique(rc_labels)))
        hold on
        boxplot(network_measures(:,imeasure,max(rc_labels(:))+1-label),agegroup_labels, 'Labels',group_names);
        ylim([min(min(network_measures(:,imeasure,:))) max(max(network_measures(:,imeasure,:)))])
        
        % remove unwanted x and y labels
        if (label+(imeasure-1)*4) ~= 1 && (label+(imeasure-1)*4) ~= 5 && (label+(imeasure-1)*4) ~= 9
            set(gca, 'YTickLabel','')
        end
        if (label+(imeasure-1)*4) < 9
            set(gca, 'XTickLabel','')
        end
        
        % adjust subplot size
        set(gca,'Position',[subplot_x(label), subplot_y(imeasure), subplot_width, subplot_height])
    end
end

%% plot comparison hitmap
all_labels = zeros(size(rc_labels,1),2*size(rc_labels,2));
all_labels(:,1:4) = rc_labels==3;

% find top N NDI regions
NDI_top_N = zeros(size(W,1),4);
all_top_N = [];
top_N = 10;
for group = 1:4
    [aa,bb] = sort(NDI_scores(:,group),'descend');
    NDI_top_N(bb(1:top_N),group) = 1;
    all_top_N = [all_top_N; bb(1:top_N)];
end
all_labels(:,5:8) = 2*NDI_top_N;

% define group names for plotting
plot_group_names = {};
for ii = 1:4
    plot_group_names(ii) = strcat('RC: ',group_names(ii));
    plot_group_names(ii+4) = strcat('NDI: ',group_names(ii));
end

% define colors and plot
colours = [[255 255 255];[86 180 233];[213 94 0]]/255.;
get_hitmap(all_labels,plot_group_names,'Labels_170_NKI_with_header.csv',colours,0)

%% NDI of RC
figure
% subplot adjustments
subplot_x = [0.1 0.4 0.7];
for rc_sub = 1:length(unique(rc_labels(:)))
   subplot(1,3,rc_sub)
   
   % gather data
   bdata = [];
   group_idx = [];
   for igroup = 1:length(group_names)
       bdata = [bdata; NDI_scores(rc_labels(:,igroup)==rc_sub,igroup)];
       group_idx= [group_idx; igroup*ones(sum(rc_labels(:,igroup)==rc_sub),1)];
   end
   
   % create plot and labels
   boxplot(bdata, group_idx, 'Labels',group_names)
   ylim([min(NDI_scores(:)), max(NDI_scores(:))]);
   title(subnetwork_names{4-rc_sub})
   if rc_sub == 1
       ylabel('NDI')
   else
       set(gca, 'YTickLabel','')
   end
   % subplot adjustments
   set(gca,'Position',[subplot_x(rc_sub), 0.1, 0.28, 0.75])
end

%% plot comparison between groups -> Spearman correlation values
fig=figure;
set(fig, 'Position', [100 100 1600 800])
hold on
all_correlations = [];
% compare group ii
for ii = 2:size(NDI_scores,2)
    % to group jj
    for jj = 1:(ii-1)
        % scatter plot
        subplot(3,3,(ii-1)+(jj-1)*3)
        scatter(NDI_scores(:,ii), NDI_scores(:,jj),'f')
        
        % calculate spearman correlation
        r = corr(NDI_scores(:,ii), NDI_scores(:,jj),'type','Spearman');
        all_correlations = [all_correlations; r];
        text(0.00000005,0.1,sprintf('r=%0.3f',r));
        
        % add labels and adjust plot
        if ((ii-1)+(jj-1)*3) == 1 || ((ii-1)+(jj-1)*3) == 5 || ((ii-1)+(jj-1)*3) == 9
            xlabel(group_names(ii))
            ylabel(group_names(jj))
        else
            set(gca,'YTickLabel',[],'XTickLabel',[]);
        end
        xlim([min(NDI_scores(NDI_scores>0)), 1])
        ylim([min(NDI_scores(NDI_scores>0)), 1])
        set(gca,'YScale','log')
        set(gca,'XScale','log')
    end
end


%% Generate confusion matrices

% size(NDI_labels) (170,4): 4 Tiers
% size(rc_labels) (170,4): 3 subnetworks

figure('Position', [10 10 1550 1000])

%%%%%%
% NDI
%%%%%%
% prep adding numbers in the imagesc plot
num_labels = length(unique(NDI_labels(:)));
x = repmat(1:num_labels, num_labels,1);
y = x';

% initialize image
confusion = zeros(num_labels);
% run through age group i
for igroup = 1:size(NDI_labels,2)
    % compare with age group j
    for jgroup = (igroup+1):size(NDI_labels,2)
        % for each label in age group i
        for ilabel = 1:num_labels
            % find the labels in age group j
            idx = NDI_labels(:,igroup)==ilabel;
            values = NDI_labels(idx,jgroup);
            % determine confusion
            for jlabel = 1:num_labels
                confusion(jlabel, ilabel) = sum(values==jlabel); 
            end
        end
        
        % plot
        subplot(size(NDI_labels,2), size(NDI_labels,2), (igroup-1)*(size(NDI_labels,2)) + (jgroup))
        % normalize confusion matrix to represent percent
        confusion = round(100*confusion./repmat(sum(confusion,1),size(confusion,2),1));
        imagesc(confusion)
        colormap(viridis)
        xlabel(group_names{igroup})
        ylabel(group_names{jgroup})
        set(gca, 'YTick', 1:4, 'XTick', 1:4)

        % add the numbers
        t = num2cell(confusion);
        t = cellfun(@num2str, t, 'UniformOutput', false);
        text(x(:), y(:), t, 'HorizontalAlignment', 'Center','Color',[1,1,1],'FontWeight','bold')
    end
end

%%%%%%
% RC
%%%%%%

% adjust order of labelling for RC framework
rc_labels = 4-rc_labels;
rc_names = {'RC','F','S'};

% prep adding numbers in the imagesc plot
num_labels = length(unique(rc_labels(:)));
x = repmat(1:num_labels, num_labels,1);
y = x';

% initialize image
confusion = zeros(num_labels);
for igroup = 1:size(rc_labels,2)
    % run through age group i
    for jgroup = (igroup+1):size(rc_labels,2)
        % compare with age group j
        for ilabel = 1:num_labels
            % find the labels in age group j
            idx = rc_labels(:,igroup)==ilabel;
            values = rc_labels(idx,jgroup);
            % determine confusion
            for jlabel = 1:num_labels
                confusion(jlabel, ilabel) = sum(values==jlabel); 
            end
        end
        
        % plot
        subplot(size(rc_labels,2), size(rc_labels,2), (jgroup-1)*(size(rc_labels,2)) + (igroup))
        % normalize confusion matrix to represent percent
        confusion = round(100*confusion./repmat(sum(confusion,1),size(confusion,2),1));
        imagesc(confusion)
        colormap(viridis)
        xlabel(group_names{igroup})
        ylabel(group_names{jgroup})

        % add the numbers
        t = num2cell(confusion);
        t = cellfun(@num2str, t, 'UniformOutput', false);
        text(x(:), y(:), t, 'HorizontalAlignment', 'Center','Color',[1,1,1],'FontWeight','bold')
        set(gca, 'YTickLabel',rc_names,'XTickLabel',rc_names)
    end
end
