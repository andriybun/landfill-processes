% Run the HYDRUS-1D model with new parameters and initial and boundary conditions and get simulated time series of soil water contents
function sims = runH1D(x,Extra)


% Modify SELECTOR.IN
ModifySelectorIn(x,Extra)

% Modify PROFILE.DAT
Extra.initial(:,3) = x(7);
ModifyProfileDat(Extra)

% Modify ATMOSPH.IN
Extra.boundcon(:,7) = x(7);
ModifyAtmosphIn(Extra)

% Run HYDRUS-1D
cd(Extra.subdir)
[status,output] = dos('H1D_CALC.EXE');
cd(Extra.workdir)

% Read OBS_NODE.OUT
sims = ReadObsNodeOut(Extra);
