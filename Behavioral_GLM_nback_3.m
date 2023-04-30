clear; clc; 
flag = 2; 
dm_matrix = 3; 

load 'schedule_testing.mat'; % schedule used in task with pertinent information. 

gains= double(baits); % Load outcomes (rewards/no rewards) from schedule


n_back = 3; % number of trials back used to predict current choices 

subs = 1:28; 
subs = subs(subs~= 7); 
subs = subs(subs~= 15);
subs = subs(subs~= 9);
subs = subs(subs~= 18);
subs = subs(subs~= 27 ); 
subs = subs(subs~= 24 ); 
subs = subs(subs~= 23 ); 

for iS=1:length(subs)
    
    load([pwd,'/choices/choices_',num2str(subs(iS)),'.mat']); % load subject choices
    
    
    c1back =toeplitz(double(choices == 1),  zeros(n_back,1)); % past history of choices
    
    r1back=toeplitz(gains(:,1),zeros(n_back,1)); % past history of Stim 1 rewards
    r2back=toeplitz(gains(:,2),zeros(n_back,1)); % past history of Stim 2 rewards
    diff_mags = mags(:,1) - mags(:,2);
    
    
    ch_out1 = gains(:,1); % chosen Stim 1 rewards
    ch_out2 = gains(:,2); % chosen Stim 2 rewards
    ch_out1(choices==2)=0;    
    ch_out2(choices==1)=0; 
    cnt_out1 = gains(:,1);
    cnt_out2 = gains(:,2); 
    cnt_out1(choices==1)=0; % counterfactual Stim 1 rewards
    cnt_out2(choices==2)=0;% counterfactual Stim 2 rewards
    
    ch_out1back = toeplitz(ch_out1, zeros(n_back,1));
    ch_out2back = toeplitz(ch_out2, zeros(n_back,1));
    cnt_out1back = toeplitz(cnt_out1, zeros(n_back,1));
    cnt_out2back = toeplitz(cnt_out2, zeros(n_back,1));
    
    if flag==1
        dm = [r1back r2back]; % this one is past rewards independent of choices
    elseif flag == 2 
        dm = [ch_out1back cnt_out1back ch_out2back cnt_out2back]; % this one has separate chosen and counterfactual outcomes for each past stim choice
    end
    
    dm=dm(1:end-1,:); % predict using outcomes from trials 1 to end-1
    
    
    if dm_matrix == 1 
    
        dm = [dm diff_mags(2:end)]; % add current trial reward magnitude difference to design matrix
        
    elseif dm_matrix == 2
        
        dm = [dm, normalize(mags(2:end,1)),normalize(mags(2:end,2))]; 
        
    elseif dm_matrix == 3 
        
        dm = [pEst,pEst2, normalize(mags(:,1)),normalize(mags(:,2))]; 
    end 

%     
    if dm_matrix == 3
        
       data1 = (choices(1:end) ==1);
       
    else
        
    data1=(choices(2:end)==1); % explain choices from trial 2 to end
%     
    end 
    [b1,dev1,stats1]=glmfit(dm,data1,'binomial','link','logit'); % Estimate logistic model 
   
    
    b(iS,:)=b1; 
end

% compute t-statistics and p-values for the group 
const=b(:,1);
b(:,1)=[];
t=mean(b)./std(b).*sqrt(length(subs)-1);
p = t_to_p(t,length(subs)-1);

set(0,'DefaultAxesTitleFontWeight','normal');

if flag == 1 
    figure('Renderer', 'painters', 'Position', [10 10 900 600])  
    imagesc(corr(dm))
    colorbar
    title 'Feedback Independent of Choice' 
    xlabel('NBack')
    xticks([1 2 3 4 5 6 ])
    xticklabels({'Face 1', 'Face 2','Face 3','House 1', 'House 2', 'House 3'})
    ylabel('NBack') 
    yticklabels({'Face 1', 'Face 2','Face 3','House 1', 'House 2', 'House 3'})
    set(gca,'FontSize',18); 
    saveas(figure(1),['DM_PC_IR_Corr'],'png')
    
    figure('Renderer', 'painters', 'Position', [10 10 900 600]) 
    boxplot(b)
    xlabel('Betas NBack')
    ylabel('Weights') 
    title 'Feedback Independent of Choice' 
    set(gca,'FontSize',18); 
    xticks([1 2 3 4 5 6 ])
    yline(0,'-.k','LineWidth',1.5);
    xticklabels({'Face 1', 'Face 2','Face 3','House 1', 'House 2', 'House 3'})
    set(gca,'linewidth',1.5)
    box off
    saveas(figure(2),['Betas_L_PC_IR'],'png')
    
elseif flag == 2 
    figure('Renderer', 'painters', 'Position', [10 10 900 600]) 
    imagesc(corr(dm))
    colorbar
    title 'Outcome and Counterfactuals Influence on Choice' 
    xticks([1 2 3 4 5 6 7 8 9 10 11 12 ])
    xticklabels({'F1', 'F2','F3', 'F1 CF', 'F2 CF','F3 CF',...
        'H1', 'H2', 'H3','H1 CF', 'H2 CF', 'H3 CF' })
    
    
    yticks([1 2 3 4 5 6 7 8 9 10 11 12 ])
    yticklabels({'F1', 'F2','F3', 'F1 CF', 'F2 CF','F3 CF',...
        'H1', 'H2', 'H3','H1 CF', 'H2 CF', 'H3 CF' })
    
    set(gca,'FontSize',20); 
    saveas(figure(1),['DM_SR_CO_Corr'],'png')
    
    %figure('Renderer', 'painters', 'Position', [10 10 900 600]) 
    figure()
    boxplot(b)
    xlabel('Nback')
    ylabel('Beta Weights') 
    title 'Trial Influence on Choice'  
    set(gca,'FontSize',18);
    set(gca,'TickDir','out'); % The only other option is 'in
    yline(0,'--k','LineWidth',2);
    set(gca,'linewidth',2)
    ylim([-5.5 5.5])
    xticks([1 2 3 4 5 6 7 8 9 10 11 12 ])
   xticklabels({'F1', 'F2','F3', 'F1 CF', 'F2 CF','F3 CF',...
        'H1', 'H2', 'H3','H1 CF', 'H2 CF', 'H3 CF' })
    xtickangle(45)
    box off 
    saveas(figure(1),['learner_regression'],'svg')
end 

    


% Corr Matrix 
if dm_matrix == 1 
    figure()
    dm_mod = dm(:,1:6); 
    imagesc(corr(dm_mod))
    colorbar
    title 'Design Matrix, Difference Mags' 
    set(gca,'FontSize',16); 
    saveas(figure(1),['DM_Diff_Mags_Corr_3'],'png')
    
elseif dm_matrix == 2 
    figure()
    dm_mod = dm(:,1:6); 
    imagesc(cov(dm_mod))
    colorbar
    title 'Design Matrix, Normalized Mags' 
    set(gca,'FontSize',16); 
    saveas(figure(1),['DM_Norm_Mags_Corr_3'],'png')
    
elseif dm_matrix == 3 
    figure('Renderer', 'painters', 'Position', [10 10 900 600]) 
    imagesc(corr(dm))
    colorbar
    title 'Probablity Estimates and Normalized Mags' 
    xticks([1 2 3 4])
    xticklabels({'pEst Face', 'pEst House','Face Mag','House Mag'})
    yticks([1 2 3 4])
    yticklabels({'pEst Face', 'pEst House','Face Mag','House Mag'})
    set(gca,'FontSize',18); 
    saveas(figure(1),['DM_Pest_Mags_Betas_Corr_3'],'png')
    
end 
    
Betas 
if dm_matrix == 1
    figure('Renderer', 'painters', 'Position', [10 10 900 600]) 
    err = std(b)/sqrt(length(b)); 
    e = errorbar(mean(b),err,'o'); 
    e.Color = 'k';
    e.MarkerSize = 10;
    e.CapSize = 15;
    ylim([-4 4])
    xlim([0 8])
    title 'Betas Original Design Matrix, Mag Diff' 
    xlabel('Betas')
    ylabel('Weights') 
    set(gca,'FontSize',18); 
    saveas(figure(2),['DM_Diff_Mags_Betas_3'],'png')
    
elseif dm_matrix == 2 
    figure()
    err = std(b)/sqrt(length(b)); 
    e = errorbar(mean(b),err,'o');
    e.Color = 'k';
    e.MarkerSize = 10;
    e.CapSize = 15;
    ylim([-4 4])
    xlim([0 9])
    xlabel('Betas')
    ylabel('Weights') 
    title 'Betas Normalized Mags' 
    set(gca,'FontSize',14); 
    saveas(figure(2),['DM_Norm_Mags_Betas_3'],'png')
    
elseif dm_matrix == 3
    figure()
    boxplot(b)
    xlabel([])
    ylabel('Beta Weights') 
    title 'pEst & Mag Choice Influence' 
    yline(0,'--k','LineWidth',2);
    set(gca,'linewidth',2)
    set(gca,'FontSize',18); 
    xticks([1 2 3 4])
    xticklabels({'pEst Face', 'pEst House','Face Mag','House Mag'})
    box off
    set(gca,'TickDir','out'); % The only other option is 'in'
    ylim([-12 12])

    saveas(figure(2),['DM_Pest_Mags_Betas_3'],'svg')
    
end 
