function plotDREAM(fn)
close all; pname = '../MSWS1'; addpath(genpath(pname)); load(fn); isubr = 4;  isubc = 4;
Pnames = {'km(hyd)' 'km(meth)' 'km(sulph)' 'km(anamm)' 'Ks(meth)' 'Ks(sulp)' 'Cxm' 'Cxs' 'Sigma_1' 'Sigma_2' 'Sigma_3' 'Sigma_4' 'Sigma_5' 'Sigma_6'};

% FIGURE 1: Plot the evolution of archive Z
Zind=max(find(Z(:,1)>0)); %the end-value of Z
figure;
for i = 1:MCMCPar.n;
    subplot(isubr,isubc,i);
    plot(Z(1:Zind,i),'k.');
    
    ylabel(Pnames{i},'FontSize',10)
    xlabel('Sample number in Z','FontSize',10)
    set(gca,'FontSize',16);
%     axis([0 Zind ParRange.minn(i) ParRange.maxn(i)])
end
%also plot the simulation with the minimum objective function/maximum probability
[minZ minZidx] = max(Z(MCMCPar.m0+1:Zind,MCMCPar.n+1));
minZidx=minZidx+MCMCPar.m0;  %the index corresponding to the parameter set with the highest probability
for i = 1:MCMCPar.n;
    subplot(isubr,isubc,i);hold on;
    plot(minZidx,Z(minZidx,i),'ro');
end

% FIGURE 2: Plot the evolution of the individual sequences (the Markov chains)
figure;
for i = 1:MCMCPar.n;
    Sind=max(find(Sequences(:,1,1)>0));
    Seq=reshape(Sequences(1:Sind,i,:),Sind,3);
    subplot(isubr,isubc,i);
    plot(Seq,'.')
    ylabel(Pnames{i},'FontSize',10)
    xlabel('Sample number in seq','FontSize',10)
    set(gca,'FontSize',16);
%     axis([0 Sind ParRange.minn(i) ParRange.maxn(i)])
end

% FIGURE 3: Plot histograms (used to calculate mean and standard deviation)
figure;
for i = 1:MCMCPar.n
    subplot(isubr,isubc,i);
    hist(Z(round(0.75.*Zind):Zind,i),10)
    xlabel(Pnames{i},'FontSize',16)
    set(gca,'FontSize',16);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[0.5 0.5 0.5])
%     xlim([ParRange.minn(i) ParRange.maxn(i)])
end

% FIGURE 4: Plot the Gelman-Rubin Convergence criterion (if criterion <1.2 for all parameters, convergence generally is declared)
figure;hold on;
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

% FIGURE 5: Plot optimum solution
SimPar_opt = Z(minZidx,1:MCMCPar.n); Z(minZidx,:)
[~, ~, Extra.ORI, ~, ~] = initialize('Pmatrix.csv');
[~, sim] = integrate_bioreactor(SimPar_opt,Extra);
plot_integration(sim.t, sim.MT, sim.t_D_o, sim.C_D_o, Extra.Meas, Extra.Pm);

% FIGURE 6: Plot optimum solution
figure;
plot(Z(1:Zind,end),'x');
end
