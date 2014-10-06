function [pv, MobileC, LongTermC, mRemaining] = Rainfall(...
    t, dt, iT, pv, mRemaining, RainInfo, MobileC, LongTermC, SpeciesInfo, ModelParams, Const)
    %

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
            error('Not implemented!');
        else
            % Increase volume of water leaving the system at future times
%             MobileC = MobileC.AddSolute(qOutAfter, ...
%                 repmat(RainInfo.concentration(iT), [1, nCalcT, SpeciesInfo.nFlush]), ...
%                 iT, iT:iTend, SpeciesInfo.iFlush);
            MobileC = MobileC.AddSolute(qOutAfter, ...
                repmat(RainInfo.concentration(iT), [1, nCalcT, SpeciesInfo.nFlush]), ...
                1, iT:iTend, SpeciesInfo.iFlush);
            LongTermC = LongTermC.AddSolute(RainInfo.intensity(iT) - sum(qOutAfter), ...
                repmat(RainInfo.concentration(iT), [1, 1, SpeciesInfo.nFlush]), ...
                :, :, SpeciesInfo.iFlush);
            mRemaining(2, iT, :) = mRemaining(2, iT, :) + ...
                RainInfo.intensity(iT) * RainInfo.concentration(iT);
        end
    end
end