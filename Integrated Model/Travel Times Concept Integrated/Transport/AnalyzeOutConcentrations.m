function AnalyzeOutConcentrations(ModelOutput, TimeParams, ModelParams, ...
    ParameterOfInterest, iSpecies, Const, PLOT_SEPARATELY, tEnd)

    specieTitles = {...
        1,  'C(hyd)';
        2,  'H(acetate)';
        3,  'H_2CO_3';
        4,  'NH_3';
        5,  'CH_4';
        6,  'H_2O';
        7,  'H+';
        8,  'SO_4-2';
        9,  'H_2S';
        10, 'C_x(ace)';
        11, 'C_x(meth)';
        12, 'C_x(sulf)';
        13, '';
        14, '';
        15, '';
        16, '';
        17, '';
        18, '';
        19, '';
        20, '';
        21, '';
        22, ''
        };

    nT = min(ModelOutput.nT, TimeParams.numIntervals);
    tShow = true(1, nT);
    
    if ~PLOT_SEPARATELY
        iSpecies = [1:8, 9:12, 17:20];
        hMax = 4;
        vMax = 4;
        figH = figure();
        set(figH, 'Position', [400, 1, 1000, 1640]);
    else
        hMax = numel(iSpecies);
        vMax = 1;
    end
    
    for iH = 1:hMax
        for iV = 1:vMax
            iP = vMax * (iV - 1) + iH;
            
            pv = [ModelParams.totalPv - ModelParams.totalPv / (1 + ModelParams.beta); ...
                ModelParams.totalPv / (1 + ModelParams.beta)];
            emissionPotential = (pv(1) * ModelOutput.cRemaining(1, 2:nT+1, iSpecies(iP)) + ...
                pv(2) * ModelOutput.cRemaining(2, 2:nT+1, iSpecies(iP))) / ModelParams.totalPv;
            leachateConcentration = ModelOutput.mOutTotal(:, tShow, iSpecies(iP)) ./ ...
                ModelOutput.qOutTotal(1, tShow);
            
            % [cMax, tMax] = findpeaks(leachateConcentration);
            % [cMin, tMin] = findpeaks(-leachateConcentration);
            % [tMax, cMax] = RemoveDrops(tMax, cMax);
            % [tMin, cMin] = RemoveDrops(tMin, cMin);
            % cMin = -cMin;
            
            if PLOT_SEPARATELY
                fH = figure();
            else
                subplot(vMax, hMax, iP);
            end
            plot(TimeParams.t(tShow), leachateConcentration, 'Color', [1, 0.5, 0.15]);
            hold on;
            %     plot(TimeParams.t(tMax), cMax, 'Color', [0.6, 0, 0], 'LineWidth', 2);
            %     plot(TimeParams.t(tMin), cMin, 'Color', [0.6, 0, 0], 'LineWidth', 2);
            plot(TimeParams.t(tShow), emissionPotential, 'Color', [0, 0.7, 0], 'LineWidth', 2, ...
                'LineStyle', '--');
            hold off;
            
            if (nargin == 8)
                xlim([0, tEnd]);
            end
            yTop = max(max(emissionPotential) * 1.2);
            if ~RealEq(yTop, 0, Const.EPSILON)
                ylim(sort([0, yTop]));
            end
            
            if PLOT_SEPARATELY
                title(sprintf('Concentration dynamics (sp. %s (#%02d), %s = %3.2f, %s = %4.3f)', ...
                    specieTitles{iSpecies(iP), 2}, ...
                    iSpecies(iP), ...
                    '\beta', ModelParams.beta, ...
                    'k_e_x', ModelParams.kExch));
                legend('leachate concentration', 'remaining concentration');
                % legend('leachate concentration', 'concentration local maxima', ...
                %     'concentration local minima', 'remaining concentration');
                xlabel('time [days]');
                ylabel('concentration [kg/m^3]');
                fileName = sprintf('fig/Cout_vs._Crem_sp_#%02d_%s.fig', ...
                    iSpecies(iP), GenerateCharacteristicSuffix(ModelParams, ParameterOfInterest));
                hgsave(fH, fileName);
            else
                title(specieTitles{iSpecies(iP), 2});
                % xlabel('t [days]');
                % ylabel('C [kg/m^3]');
            end
        end
    end
    
    if ~PLOT_SEPARATELY
        % suptitle(sprintf('Concentration dynamics (%s = %3.2f, %s = %4.3f)', ...
        %     '\beta', ModelParams.beta, ...
        %     'k_e_x', ModelParams.kExch));
    end
    
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