function plotblim(model0,filename,options)
% PLOTBLIM Plots summaries and disgnostic plot for the BLIM
% Return a series of plots
%
% PLOTBLIM(MODEL0,FILENAME,OPTIONS)
%
% where:
% - model0 is a the outcome of a function for fitting the blim
% - filename is a optional parameters. If filename is a string the
%            plot will be save in a file with that name
% 
% PLOTBLIM(...,OPTIONS)
%
%   options.parameters: tolerance for parameter change 
%
% Authors:  Andrea Brancaccio (andrea.brancaccio@unipd.it)
if nargin<3
    options.parameters='error';
    options.residual='Brand-Altman';
end

if(size(model0.states,1)<50)
    figure(1)
    plot_knowledge_structure(model0.states)
    if nargin ==2
        exportgraphics(gca,filename,"Append",true)
    end
end
figure(2)
t1=plotparameters(model0,options.parameters);
if nargin ==2
    exportgraphics(t1,filename,"Append",true)
end
% figure(3)
% plotout(model0,options.residual)
% if nargin ==2
%     exportgraphics(gca,filename,"Append",true)
% end
figure(3)
t2=plotboot(model0);
if nargin ==2
    exportgraphics(t2,filename,"Append",true)
end
end

function t=plotparameters(model,type)
t = tiledlayout(1,2);
outliers=isoutlier(model.eta_bootstrap,'quartiles');
nexttile
if strcmp(type,"box")
    outliers=isoutlier(model.beta_bootstrap,'quartiles');
    %boxchart(model.beta_bootstrap(1-outliers)')
    beta=[];
    items=[];
    for i=1:size(outliers,1)
        beta = [beta,model.beta_bootstrap(i,outliers(1,:)==0)];
        items = [items , repmat(i,1,sum(outliers(1,:)==0))];
    end
    boxchart(items, beta)
    ylim([0,1])
    xlim([0,size(model.classes,1)+1]) %% aggiunto!!
    title('$\beta$ estimates','interpreter','latex','fontsize',14)
    xlabel('Item','interpreter','latex','fontsize',14)
    ylabel('Probability','interpreter','latex','fontsize',14)
    hold on
    plot([0,size(model.classes,1)+1],[.5,.5],'k--')
    hold off
    set(gca,'fontsize',14)
    
    nexttile
    
    eta=[];
    items=[];
    for i=1:size(outliers,1)
        eta = [eta,model.eta_bootstrap(i,outliers(1,:)==0)];
        items = [items , repmat(i,1,sum(outliers(1,:)==0))];
    end
    boxchart(items, eta)
    ylim([0,1])
    xlim([0,size(model.classes,1)+1]) %% aggiunto!!
    title('$\eta$ estimates','interpreter','latex','fontsize',14)
    xlabel('Item','interpreter','latex','fontsize',14)
    ylabel('Probability','interpreter','latex','fontsize',14)
    hold on
    plot([0,size(model.classes,1)+1],[.5,.5],'k--')
    hold off
    set(gca,'fontsize',14)
    
    
elseif strcmp(type,"error")
    
    errorbar([1:size(model.beta,1)]',model.beta,model.beta_sd,'k*' )
    ylim([0,1])
    xlim([0,size(model.classes,1)+1]) %% aggiunto!!
    title('$\beta$ estimates','interpreter','latex','fontsize',14)
    xlabel('Item','interpreter','latex','fontsize',14)
    ylabel('Probability','interpreter','latex','fontsize',14)
    hold on
    plot([0,size(model.classes,1)+1],[.5,.5],'k--')
    hold off
    set(gca,'fontsize',14)
    
    nexttile
    
    errorbar([1:size(model.eta,1)]',model.eta,model.eta_sd,'k*' )
    ylim([0,1])
    xlim([0,size(model.classes,1)+1]) %% aggiunto!!
    title('$\eta$ estimates','interpreter','latex','fontsize',14)
    xlabel('Item','interpreter','latex','fontsize',14)
    ylabel('Probability','interpreter','latex','fontsize',14)
    hold on
    plot([0,size(model.classes,1)+1],[.5,.5],'k--')
    hold off
    set(gca,'fontsize',14)
end
end

function t = plotboot(model)
if isnan(model.LR0)==true
       nexttile
    plot(sort(model.chi_bootstrap),1-(1:model.nrep)/model.nrep, 'k')
    title(sprintf('Chi Squared bootstrap distribution\n %i iteration',model.nrep),'interpreter','latex','fontsize',14)
    xlabel('bootstraped $\chi^2$','interpreter','latex','fontsize',14)
    ylabel('$p$-value','interpreter','latex','fontsize',14)
    hold on
    plot([model.chi model.chi],[0 1-sum(model.chi_bootstrap<model.chi)/model.nrep], 'k--')
    plot([0 model.chi],[model.pvalue_chi model.pvalue_chi], 'k--')
    hold off
else
    t = tiledlayout(2,1);
    nexttile
    plot(sort(model.LR_bootstrap),1-(1:model.nrep)/model.nrep, 'k')
 title(sprintf('Likelihood ratio bootstraped distribution\n %i iteration',model.nrep),'interpreter','latex','fontsize',14)
 xlabel('bootstraped LR','interpreter','latex','fontsize',14)
    ylabel('$p$-value','interpreter','latex','fontsize',14)
    hold on
    plot([model.LR0 model.LR0],[0 1-sum(model.LR_bootstrap<model.LR0)/model.nrep], 'k--')
    plot([0 model.LR0],[model.pvalue_LR model.pvalue_LR], 'k--')
    hold off
    
    nexttile
    plot(sort(model.chi_bootstrap),1-(1:model.nrep)/model.nrep, 'k')
    title(sprintf('Chi Squared bootstrap distribution\n %i iteration',model.nrep),'interpreter','latex','fontsize',14)
    xlabel('bootstraped $\chi^2$','interpreter','latex','fontsize',14)
    ylabel('$p$-value','interpreter','latex','fontsize',14)
    hold on
    plot([model.chi model.chi],[0 1-sum(model.chi_bootstrap<model.chi)/model.nrep], 'k--')
    plot([0 model.chi],[model.pvalue_chi model.pvalue_chi], 'k--')
    hold off
    
end

end


function  plotout(model,type)
if strcmp(type,"Brand-Altman",'interpreter','latex','fontsize',14)
    residuals = model.freq-model.expfreq;
    average = (model.freq+model.expfreq)/2;
    plot( average, residuals,'k*')
    title(sprintf('Brand-Altman plot'),'interpreter','latex','fontsize',14)
    ylabel('Residuals')
    xlabel('(Expected+Observed)/2','interpreter','latex','fontsize',14)
    hold on
    plot([min(average),max(average)+1],[mean(residuals),mean(residuals)],'k--')
    plot([min(average),max(average)+1],[mean(residuals)+2*std(residuals,[]),mean(residuals)+2*std(residuals,[])],'k-.')
    plot([min(average),max(average)+1],[mean(residuals)-2*std(residuals,[]),mean(residuals)-2*std(residuals,[])],'k-.')
    %loglog([0:.0001:1],[0:.0001:1])
    hold off
elseif strcmp(type,"observed vs expected")
    plot(model.expfreq./model.sz,model.freq/model.sz,'*')
    title(sprintf('Observed VS Expected'),'interpreter','latex','fontsize',14)
    xlabel('Expected proportion','interpreter','latex','fontsize',14)
    ylabel('Observed proportion','interpreter','latex','fontsize',14)
    hold on
    %loglog([0:.0001:1],[0:.0001:1])
    hold off
end

end
