function [pv, PartInfo, mRemaining] = Rainfall(...
    t, dt, iT, pv, mRemaining, RainInfo, PartInfo, SpeciesInfo, ModelParams, Const)
    %

    nPartGeneral = 2;
    iImmobile = 1;
    iMobileIni = 2;
    
    if RealGt(RainInfo.intensity(iT), 0, Const.EPSILON)
        nT = numel(t);
        tOffset = t(iT);
        tAfter = t(iT:nT) - tOffset;
        tBounds = ModelParams.LogNorm.bounds;
        % Select only those time steps that are affected by current injection
        iCalcT = (tAfter <= tBounds(2));
        nCalcT = sum(iCalcT);
        iTend = iT + nCalcT - 1;
        tAfter = tAfter(iCalcT);
    
        % Add (fresh) rainwater to the system. Change volume of liquid and revise mass and
        % concentration of solute in mobile phase
        pvMobUpd = pv(2) + RainInfo.intensity(iT);
        pv(2) = pvMobUpd;
        % Every input impulse of water will cause (log-normal) response at the outlet. All the
        % outflow during a given time step is considered as a particle with unique travel time
        if iT == nT
            lnPdf = ModelParams.LogNorm.PdfDelayed([tAfter, tAfter+dt]);
            lnPdf = lnPdf(1);
        else
            lnPdf = ModelParams.LogNorm.PdfDelayed(tAfter);
        end
        qOutAfter = RainInfo.intensity(iT) * lnPdf;
        if (ModelParams.DO_RECIRCULATION)
            if ((iT > 1) && (~RealEq(PartInfo.GetVolume(1, iT-1), 0, Const.EPSILON)))
                cOutPrevStep = mOutTotal(1, iT-1, SpeciesInfo.iFlush) / PartInfo.GetVolume(1, iT-1);
                mRemaining(2, iT, SpeciesInfo.iFlush) = mRemaining(2, iT, SpeciesInfo.iFlush) + ...
                    RainInfo.intensity(iT) * cOutPrevStep;
                PartInfo = PartInfo.AddSolute([qOutAfter, RainInfo.intensity(iT) - sum(qOutAfter)], ...
                    repmat(cOutPrevStep, [1, nCalcT+1, 1]), ...
                    :, [iMobileIni, nPartGeneral + (iT:iTend)], SpeciesInfo.iFlush);
            end
        else
            % Increase volume of water leaving the system at future times
            PartInfo = PartInfo.AddSolute([qOutAfter, RainInfo.intensity(iT) - sum(qOutAfter)], ...
                repmat(RainInfo.concentration(iT), [1, nCalcT+1, SpeciesInfo.nFlush]), ...
                :, [iMobileIni, nPartGeneral + (iT:iTend)], SpeciesInfo.iFlush);
            mRemaining(2, iT, :) = mRemaining(2, iT, :) + ...
                RainInfo.intensity(iT) * RainInfo.concentration(iT);
        end
    end
end