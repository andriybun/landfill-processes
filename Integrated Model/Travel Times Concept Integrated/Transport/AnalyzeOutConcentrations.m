function AnalyzeOutConcentrations(ModelOutput, TimeParams, ModelParams, ...
    ParameterOfInterest, iSpecies)

    nT = min(ModelOutput.nT, TimeParams.numIntervals);
    tShow = true(1, nT);
    
    pv = [ModelParams.totalPv - ModelParams.totalPv / (1 + ModelParams.beta); ...
        ModelParams.totalPv / (1 + ModelParams.beta)];
    emissionPotential = (pv(1) * ModelOutput.cRemaining(1, 2:nT+1, iSpecies) + ...
        pv(2) * ModelOutput.cRemaining(2, 2:nT+1, iSpecies)) / ModelParams.totalPv;
    leachateConcentration = ModelOutput.mOutTotal(:, tShow, iSpecies) ./ ...
        ModelOutput.qOutTotal(1, tShow);
    
    [cMax, tMax] = findpeaks(leachateConcentration);
    [cMin, tMin] = findpeaks(-leachateConcentration);
    
    [tMax, cMax] = RemoveDrops(tMax, cMax);
    [tMin, cMin] = RemoveDrops(tMin, cMin);
    
    cMin = -cMin;
    
    fH = figure();
    plot(TimeParams.t(tShow), leachateConcentration, 'Color', [1, 0.5, 0.15]);
    hold on;
    plot(TimeParams.t(tMax), cMax, 'Color', [0.6, 0, 0], 'LineWidth', 2);
    plot(TimeParams.t(tMin), cMin, 'Color', [0.6, 0, 0], 'LineWidth', 2);
    plot(TimeParams.t(tShow), emissionPotential, 'Color', [0, 0.7, 0], 'LineWidth', 2);
    hold off;
    xlim([0, 200]);
    ylim([0, max(max(emissionPotential) * 1.2)]);
    
    title(sprintf('Concentration dynamics (sp. #%02d, %s = %3.2f, %s = %4.3f)', ...
        iSpecies, ...
        '\beta', ModelParams.beta, ...
        'k_e_x', ModelParams.kExch));    % 'k_e_x_P_a_r_t', ModelParams.kExchPart)
    legend('leachate concentration', 'concentration local maxima', ...
        'concentration local minima', 'remaining concentration');
    xlabel('time [days]');
    ylabel('concentration [kg/m^3]');
    
    fileName = sprintf('fig/Cout_vs._Crem_%s_sp_#%02d.fig', ...
        GenerateCharacteristicSuffix(ModelParams, ParameterOfInterest), iSpecies);
    hgsave(fH, fileName);
    
    return
    
    function [xIn, yIn] = RemoveDrops(xIn, yIn)
        nEl = numel(xIn);
        i = 2;
        while (i < nEl - 1)
            if ((yIn(i) < yIn(i-1)) && (yIn(i) < yIn(i+1)))
                xIn(i) = [];
                yIn(i) = [];
                nEl = nEl - 1;
            else
                i = i + 1;
            end
        end
    end
end