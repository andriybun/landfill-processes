function plotDREAM(MCMCPar,Z,ParRange,Sequences,output,ExtraPar,Measurement)
close all

ExtraPar.ORI = initialize_ORI(ExtraPar.Comp,0); tm = ExtraPar.tm;
Pnames = ExtraPar.Pnames; Dnames = ExtraPar.Dnames; MeasData = Measurement.MeasData;

%% And plot some results
isubr = 4;  isubc = 4;

%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% FIGURE 1: Plot the evolution of archive Z

Zind=max(find(Z(:,1)>0)); %the end-value of Z
figure(1);
for i = 1:MCMCPar.n;
    subplot(isubr,isubc,i);
    plot(Z(1:Zind,i),'k.');
    
    ylabel(Pnames{i},'FontSize',10)
    xlabel('Sample number in Z','FontSize',10)
    set(gca,'FontSize',16);
    axis([0 Zind ParRange.minn(i) ParRange.maxn(i)])
end

%also plot the simulation with the minimum objective function/maximum probability
[minZ minZidx] = max(Z(MCMCPar.m0+1:Zind,MCMCPar.n+1));
minZidx=minZidx+MCMCPar.m0;  %the index corresponding to the parameter set with the highest probability
for i = 1:MCMCPar.n;
    subplot(isubr,isubc,i);hold on;
    plot(minZidx,Z(minZidx,i),'ro');
end
%---------------------------------------------------------------------------------------------------

%--------------------------------------------------------------------------------------------------- 
% FIGURE 2: Plot the evolution of the individual sequences (the Markov chains)
 
figure(2);
for i = 1:MCMCPar.n;
    Sind=max(find(Sequences(:,1,1)>0));
    Seq=reshape(Sequences(1:Sind,i,:),Sind,3);
    
    subplot(isubr,isubc,i);
    plot(Seq,'.')
    
    ylabel(Pnames{i},'FontSize',10)
    xlabel('Sample number in seq','FontSize',10)
    set(gca,'FontSize',16);
    axis([0 Sind ParRange.minn(i) ParRange.maxn(i)])
end
%---------------------------------------------------------------------------------------------------

%--------------------------------------------------------------------------------------------------- 
% FIGURE 3: Plot histograms (used to calculate mean and standard deviation)

figure(3);
for i = 1:MCMCPar.n
    subplot(isubr,isubc,i);
    hist(Z(round(0.80.*Zind):Zind,i),10)
    
    xlabel(Pnames{i},'FontSize',16)
    set(gca,'FontSize',16);
%     xlim([ParRange.minn(i) ParRange.maxn(i)])
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[0.5 0.5 0.5])
end
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% FIGURE 4: Plot the Gelman-Rubin Convergence criterion
%           (if criterion <1.2 for all parameters, convergence generally is declared)

figure(4);hold on;
GR_ind = max(find(output.R_stat(:,1)>0));

plot(output.R_stat(1:GR_ind,1),output.R_stat(1:GR_ind,2:end),'.-');

x = 1:(output.R_stat(GR_ind)+200);
y = ones(size(x)).*1.2;

plot(x,y,'k--');

axis([-Inf Inf 0 10])
hl=findobj(gca,'Type','Line');
legend(hl(length(Pnames)+1:2),Pnames)
xlabel('number of model evaluations')
ylabel('Gelman-Rubin convergence criterion')
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% FIGURE 5: Plot optimum solution

SimPar_opt = Z(minZidx,1:MCMCPar.n)
SimData = integrate_bioreactor(SimPar_opt,ExtraPar);
N = size(Measurement.MeasData,1);
Residual = Measurement.MeasData-SimData;
RMSE = sqrt(sum(Residual.^2)./(N-MCMCPar.n))

figure(5);
for i = 1:length(ExtraPar.tm)
    subplot(2,3,i)
    plot(tm{i},MeasData(1:length(tm{i})),'x','MarkerSize',7); hold on;
    plot(tm{i},SimData(1:length(tm{i})),'r-','LineWidth',5,'MarkerSize',10)
    MeasData(1:length(tm{i})) = []; SimData(1:length(tm{i})) = [];
    
    xlabel('Time (days)','FontSize',13)
%     ylabel('CO2+CH4 (mol)','FontSize',13)
    title(Dnames{i},'FontSize',15);
    hl=findobj(gca,'Type','Line');
    set(gca,'FontSize',15)
    legend([hl(2) hl(1)],{'Experiment','Simulation'},'Location','SouthEast','FontSize',10)
end

end
