% compute average decoding accuracy and perform cluster-based permutation

% analysis.

% Note: Randomization routine included in this analysis (e.g., random integer generator)

% can produce stats slightly different from those reported in the paper.

% Carlos Carrasco

sublist = 1:28; %num_subjects 
sublist = sublist(sublist~= 7 ); 
sublist = sublist(sublist~= 15 ); 
sublist = sublist(sublist~= 9 ); 
sublist = sublist(sublist~= 18 ); 
sublist = sublist(sublist~= 27 ); 
sublist = sublist(sublist~= 24 ); 
sublist = sublist(sublist~= 23 ); 

addpath([pwd,'/Alpha_Choice_Results/']);
set(0,'DefaultAxesTitleFontWeight','normal'); % thin title weight 

Nsub = length(sublist);

Nblock = 3; % cross-validation
Nitr = 100; % iteration
Ntp = 325; % # of time points
NBins = 2; % # of location bin

tm = -500:20:5980;

%create empty matrix
AverageAccuracy = nan(Nsub,Ntp);
    
for sub = 1:Nsub
    DecodingAccuracy = nan(Ntp,Nblock,Nitr);
    
     %% load SVM_ECOC output files
    fileLocation = [pwd, '/Alpha_Choice_Results/'];
    
    readThis =strcat(fileLocation,'CA_FvH_Choice_Alpha_Decoding_',num2str(sublist(sub)),'.mat');

    load(readThis)
     

    % prediciton from SVM-ECOC model
    svmPrediction = squeeze(svmECOC.modelPredict);
    tstTargets = squeeze(svmECOC.targets);
    clear svmECOC
    
    % compute decoding accuracy of each decoding trial
    for block = 1:Nblock
        for itr = 1:Nitr
            for tp = 1:Ntp  

                prediction = squeeze(svmPrediction(itr,tp,block,:)); % this is predictions from models
                TrueAnswer = squeeze(tstTargets(itr,tp,block,:)); % this is predictions from models
                Err = TrueAnswer - prediction;
                ACC = mean(Err==0);
                DecodingAccuracy(tp,block,itr) = ACC; % average decoding accuracy

            end
        end
    end
      
     %average across block and iterations
     grandAvg = squeeze(mean(mean(DecodingAccuracy,2),3));
    

     AverageAccuracy(sub,:) =grandAvg; % putting sub data in matrix 
     
end %End of subject

%compute average accuracy across participants and SE of the mean.
subAverage = squeeze(mean(AverageAccuracy,1)); 
seAverage = squeeze(std(AverageAccuracy,1))/sqrt(Nsub); 

%% do cluster mass analyses
 releventTime = 26:325; % from 220 ms - 1496 ms

% t-test at each relevent time point
Ps = nan(2,length(releventTime));
    for i = 1:length(releventTime)
        tp = releventTime(i);
 
        [H,P,CI,STATS] =  ttest(AverageAccuracy(:,tp), 0.5, 'tail', 'right'); % one sample t-test 

        Ps(1,i) = STATS.tstat;
        Ps(2,i) = P;
    end
    
% find significant time points
candid = Ps(2,:) <= .05;

%remove orphan time points
candid_woOrphan = candid;
candid_woOrphan(1,1) = candid(1,1);
for i = 2:(length(releventTime)-1)
    
    if candid(1,i-1) == 0 && candid(1,i) ==1 && candid(1,i+1) ==0
    candid_woOrphan(1,i) = 0; 
    else
    candid_woOrphan(1,i) = candid(1,i);     
    end
    
end

% combine whole time range with relevent time & significant information
clusters = zeros(length(tm),1);
clusterT = zeros(length(tm),1);
clusters(releventTime,1) = candid_woOrphan;
clusterT(releventTime,1) = Ps(1,:);
clusterTsum = sum(Ps(1,logical(candid_woOrphan)));

%%find how many clusters are there, and compute summed T of each cluster
tmp = zeros(10,300);
cl = 0;
member = 0;
for i = 2:length(clusters)-1

        if clusters(i-1) ==0 && clusters(i) == 1 && clusters(i+1) == 1 
        cl = cl+1;
        member = member +1;
        tmp(cl,member) = i;    
        
        elseif clusters(i-1) ==1 && clusters(i) == 1 && clusters(i+1) == 0 
        member = member +1;  
        tmp(cl,member) = i;    
        member = 0;  
        elseif clusters(i-1) ==1 && clusters(i) == 1 && clusters(i+1) == 1             
        member = member +1;  
        tmp(cl,member) = i;    
        
        else
        
        end
end


HowManyClusters = cl;
a = tmp(1:cl,:); % subset significant clusters
eachCluster = a(:,logical(sum(a,1)~=0)); % cut the size at the maximum cluster 

%now, compute summed T of each cluster 
dat_clusterSumT = nan(HowManyClusters,1);
for c = 1:HowManyClusters
   dat_clusterSumT(c,1) = sum(clusterT(eachCluster(c,eachCluster(c,:) ~=0)));
end

%% do monte-carlo simulation
NPermutations = 100;

%% note: simulation takes very long time. 

load(['Permutation_CA_FvH_Choice_Decoding_onetailed_Alpha.mat'])
   

%% find critical t-value
cutOff = NPermutations - NPermutations * 0.05; %one tailed
critT = permutedT(cutOff); % t-mass of top 95% 
sigCluster = dat_clusterSumT > critT;


%% plot significant clusters 
figure()
cl=colormap(parula(50));

% draw average accuracy
accEst = squeeze(subAverage);  
% draw clusters
draw = eachCluster(sigCluster,:);
draw = sort(reshape(draw,1,size(draw,1)*size(draw,2)));
draw = draw(draw>0);

w = zeros(Ntp,1);
w(draw)=1;
a = area(1:length(tm), accEst.*w');
a.EdgeColor = 'none';
a.FaceColor = [0.8,0.8,0.8];
child = get(a,'Children');
set(child,'FaceAlpha',0.9)
% draw mean and SE of average decoding accuracy
hold on
mEI = boundedline(1:length(tm),subAverage,seAverage, 'cmap',cl(1,:),'alpha','transparency',0.35);
xlabel('Time (ms)');ylabel('Accuracy')
set(gca,'TickDir','out'); % The only other option is 'in'

ax = gca;
ax.YLim = [0.4, .7];
ax.XTick = [0 26 51 76 101 126 151 176 201 226 251 276 301 325]; 
ax.XTickLabel = {'-500','0','500','1000','1500','2000','2500','3000', ...
    '3500','4000','4500','5000','5500','6000'};
ax.XLim = [0 325]; 
yline(0.5, '--k','LineWidth',2); 
xline(26, '--k','LineWidth',2);
set(gca,'FontSize',16); 
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'linewidth',2)
set(gca,'FontWeight','normal')

box off 
  
title('Choice Decoding') 


saveas(figure(1),'Alpha_FvH_Choice','svg')
 
    
  