
% Plotting parameters
plotLineThickness= 1.1;
plotW = 1000; % width in pixels
plotH = 300; % height in pixels

lambda1range = eta1range/sqrt(T); 


%% Separate plot

for par=1:numParams
    if plotQuietly ~= 1
        % Plot
        p = figure('Renderer', 'painters', 'Position', [50 50 plotW plotH]);
    else
        p = figure('visible','off', 'Renderer', 'painters', 'Position', [50 50 plotW plotH]);% Do not plot
    end

    for nIter=1:numN
       
       subplot(1, numN, nIter)
       v = get(gca,'Position');
       set(gca,'Position',[v(1) v(2)*1.5 v(3:4)])
    
       indRef = mseIndividual(par, :);

       plot(lambda1range, msePlugInFixed(par ,: ,nIter)./indRef,'-o', 'LineWidth', plotLineThickness)
       hold on
       plot(lambda1range, msePlugInLargeRandom1(par ,: ,nIter)./indRef,'-+', 'LineWidth', plotLineThickness)
       plot(lambda1range, msePlugInLargeRandom2(par ,: ,nIter)./indRef,'-x', 'LineWidth', plotLineThickness)
       
       plot(lambda1range, mseAIC(par ,: ,nIter)./indRef,'--', 'LineWidth', plotLineThickness)
       plot(lambda1range, mseMMA(par ,: ,nIter)./indRef, ':','LineWidth', plotLineThickness)
       plot(lambda1range, mseMG(par ,: ,nIter)./indRef,'-.', 'LineWidth', plotLineThickness)
       ylim([0.71, 1.53])
       yline(1, 'LineWidth', 0.1,'color',[.94 .94 .94])
       if nIter==1
         legend('minMSE: Fixed-N', "minMSE: $\bar{N}$="+num2str(i0_1-1), ...
             "minMSE: $\bar{N}$="+num2str(i0_2-1), 'AIC', 'MMA','Mean group', 'Location', 'best', 'Interpreter', 'Latex')  
       end
       title("N="+ Nvalues(nIter))
       xlabel('\lambda_{1}')
       
    end
    a = 'Averaging estimator, ' + des{par} + ', ratio of MSE to individual estimator';
    suptitle(a)
    
    
    figureSavingName = "Figures/relativeMSE"+num2str(numReplications)+ "etaRange"+...
        num2str(min(eta1range))+ "-" + num2str(ceil(max(eta1range)))+ "N"+...
        num2str(Nvalues(nIter))+"T"+num2str(T)+desName{par} ;
    
    saveas(p,figureSavingName,'epsc');

end


%% Plot 1
figure 

for par=1:numParams
    for nIter=1:numN
       subplot(numParams, numN, nIter+numN*(par-1))
       indRef = mseIndividual(par, :);

       plot(eta1range, msePlugInFixed(par ,: ,nIter)./indRef)
       hold on
       plot(eta1range, msePlugInLargeRandom1(par ,: ,nIter)./indRef)
       plot(eta1range, msePlugInLargeRandom2(par ,: ,nIter)./indRef)
       
       plot(eta1range, mseAIC(par ,: ,nIter)./indRef) 
       plot(eta1range, mseMMA(par ,: ,nIter)./indRef) 
       plot(eta1range, mseMG(par ,: ,nIter)./indRef) 
       ylim([0.5, 1.3])
       
       legend('P: fixed', 'P: large 1', 'P: large2', 'AIC','MMA','MG', 'Location', 'northeast')
       
       if nIter==2
           title({des{par}, "N="+ Nvalues(nIter)})
       else
           title("N="+ Nvalues(nIter))
       end
    end
end
