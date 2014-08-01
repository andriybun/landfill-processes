function [SimRR] = hymodFORTRAN(Pars,Extra); 
% Runs the FORTRAN HYMOD model and returns the simulated discharge

% Go to the subdirectory of the HYMOD model
cd(Extra.subdir);
% Write the parameter values to a file Param.in
dlmwrite('Param.in',Pars,'delimiter',' ');
% Execute the model -- this model reads the current parameter values from Param.in
dos('HYMODsilent.exe');
% Load the output of the model 
SimRR=load('Q.out');
% Return to the main directory 
cd(Extra.workdir);