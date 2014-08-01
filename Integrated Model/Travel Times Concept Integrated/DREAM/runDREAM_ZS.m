% ------------- DREAM with sampling from past and snooker updates: DREAM_ZS --------------------%
%                                                                                               %
% The code presented herein is a Markov Chain Monte Carlo algorithm that runs multiple chains   %
% in parallel for efficient posterior exploration. The algorithm, entitled DREAM_(ZS) is        %
% based on the original DREAM sampling scheme, but uses sampling from an archive of past        %
% states to generate candidate points in each individual chain. Theoy and numerical examples of %
% DREAM_(ZS) have been presented in Vrugt et al. (2009). Details can also be found in           %
% Ter Braak and Vrugt (2008)                                                                    %
%                                                                                               %
% Sampling from past has three main advantages:                                                 %
% (1) Circumvents the requirement of using N = d for posterior exploration. This will speed-up  %
% convergence to a limiting distribution, especially for high-dimensional problems (large d).   %
% (2) Outlier chains do not need explicit consideration. By sampling historical states,         %
% aberrant trajectories an jump directly to the modal region at any time during the             %
% simulation. The N path ways simulated with DREAM_(ZS) therefore maintain detailed balance at  %
% every singe step in the chain.                                                                %
% (3) The transition kernel defining the jumps in each of the chains does not require           %
% information about the current states of the chains. This is of great advantage in a           %
% multi-processor environment where the N candidate points can be generated simultaneously so   %
% that each chain can evolve most efficiently on a different computer. Details of this will be  %
% given in a later publication, which should be ready within the next few months.               %
%                                                                                               %
% DREAM_(ZS) also contains a snooker updater to maximize the diversity of candidate points      %
% and generate jumps beyond parallel direction updates. Finally, DREAM_(ZS) contains subspace   %
% learning in a similar way as DREAM, to maximize the squared jumping distance between two      %
% subsequent points in each chain. This idea has been presented in Vrugt et al. (2008) and      %
% shown to significantly increase the efficiency of posterior exploration. All these options    %
% can be activated from the input file.                                                         %
%                                                                                               %
% DREAM_(ZS) developed by Jasper A. Vrugt and Cajo ter Braak                                    %
%                                                                                               %
% This algorithm has been described in:                                                         %
%                                                                                               %
%   C.J.F. ter Braak, and J.A. Vrugt, Differential Evolution Markov Chain with snooker updater  %
%       and fewer chains, Statistics and Computing, 10.1007/s11222-008-9104-9, 2008             %
%                                                                                               %
%   Vrugt, J.A., and C.J.F. ter Braak, DiffeRential Evolution Adaptive Metropolis with Sampling %
%       from the Past and Subspace Updating, SIAM journal on Optimization                       %
%                                                                                               %
%   Vrugt, J.A., and C.J.F. ter Braak, Multiple Try DiffeRential Evolution Adaptive Metropolis  %
%       for High Performance Computing, SIAM Journal on Distributed Computing                   %
%                                                                                               %
% For more information please read:                                                             %
%                                                                                               %
%   Vrugt J.A., H.V. Gupta, W. Bouten and S. Sorooshian, A Shuffled Complex Evolution           %
%       Metropolis algorithm for optimization and uncertainty assessment of hydrologic model    %
%       parameters, Water Resour. Res., 39 (8), 1201, doi:10.1029/2002WR001642, 2003.           %
%                                                                                               %
%   ter Braak, C.J.F., A Markov Chain Monte Carlo version of the genetic algorithm Differential %
%       Evolution: easy Bayesian computing for real parameter spaces, Stat. Comput., 16,        %
%       239 - 249, doi:10.1007/s11222-006-8769-1, 2006.                                         %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson, Treatment of      %
%       input uncertainty in hydrologic modeling: Doing hydrology backward using Markov         %
%       chain Monte Carlo, Water Resour. Res., 44, W00B09, doi:10.1029/2007WR006720, 2008.      %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,       %
%       Accelerating Markov chain Monte Carlo simulation by self adaptive differential          %
%       evolution with randomized subspace sampling, International Journal of Nonlinear         %
%       Sciences and Numerical Simulation, In Press. 2009.                                      %
%                                                                                               %
% Copyright (c) 2008, Los Alamos National Security, LLC                                         %
%                                                                                               %
% All rights reserved.                                                                          %
%                                                                                               %
% Copyright 2008. Los Alamos National Security, LLC. This software was produced under U.S.      %
% Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is     %
% operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S.     %
% Government has rights to use, reproduce, and distribute this software.                        %
%                                                                                               %
% NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES A NY WARRANTY, EXPRESS OR  %
% IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE. If software is modified to    %
% produce derivative works, such modified software should be clearly marked, so as not to       %
% confuse it with the version available from LANL.                                              %
%                                                                                               %
% Additionally, redistribution and use in source and binary forms, with or without              %
% modification, are permitted provided that the following conditions are met:                   %
% � Redistributions of source code must retain the above copyright notice, this list of         %
%   conditions and the following disclaimer.                                                    %
% � Redistributions in binary form must reproduce the above copyright notice, this list of      %
%   conditions and the following disclaimer in the documentation and/or other materials         %
%   provided with the distribution.                                                             %
% � Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL %
%   the U.S. Government, nor the names of its contributors may be used to endorse or promote    %
%   products derived from this software without specific prior written permission.              %
%                                                                                               %
% THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND   %
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES      %
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS %
% ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, %
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF   %
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)        %
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT %
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,       %
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                            %
%                                                                                               %
%                                                                                               %
% Copyright (c) 2008, Los Alamos National Security, LLC                                         %
%                                                                                               %
% Written by Jasper A. Vrugt: vrugt@lanl.gov                                                    %
%                                                                                               %
% Version 0.5: January 2009                                                                     %
% Version 1.0: April 2011         Maintenance update, explicit treatment of prior distribution  %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %

% Different test examples from SIAM paper
% example 1: n-dimensional Gaussian distribution
% example 2: multivariate student t distribution
% example 3: n-dimensional banana shaped Gaussian distribution
% example 4: n-dimensional multimodal mixture distribution
% example 5: real-world example using hymod rainfall - runoff model (HYMOD code in MATLAB)
% example 6: real-world example using hymod rainfall - runoff model (HYMOD code in FORTRAN)
% example 7: rainfall-runoff model with generalized log-likelihood function
% example 8: HYDRUS-1D soil hydraulic model: using prior information on soil hydraulic parameters


example = 4;

if example == 1,    % n-dimensional Gaussian distribution

    % Problem specific parameter settings
    MCMCPar.n = 100;                        % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 300000;                 % Maximum number of function evaluations
    MCMCPar.parallelUpdate = 1.0;           % Fraction of parallel direction updates

    % Recommended parameter settings
    MCMCPar.seq = 3;                        % Number of Markov Chains / sequences
    MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.nCR = 3;                        % Number of crossover values used
    MCMCPar.m0 = 10 * MCMCPar.n;            % Initial size of Z
    MCMCPar.k = 10;                         % Thinning parameter for appending X to Z
    MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
    MCMCPar.steps = 500;                    % Number of steps before calculating convergence diagnostics
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % --------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of probability of CR values
    % --------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage -------------------
    Extra.reduced_sample_collection = 'Yes';% Only collect thinned samples
    Extra.T = 1e4;                          % Every Tth sample is collected
    % --------------------------------------------------------------------------------------------

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [-5 * ones(1,MCMCPar.n)]; ParRange.maxn = [15 * ones(1,MCMCPar.n)];
    %ParRange.minn = [9.9 * ones(1,MCMCPar.n)]; ParRange.maxn = [10 * ones(1,MCMCPar.n)];
    % Define the boundary handling
    Extra.BoundHandling = 'None';

    % ---------------------- Define covariance matrix ---------------------
    % Construct the dxd correlation matrix
    A = 0.5*eye(MCMCPar.n) + 0.5*ones(MCMCPar.n);
    % Rescale to variance-covariance matrix of interest
    for i=1:MCMCPar.n
        for j=1:MCMCPar.n
            C(i,j) = A(i,j)*sqrt(i*j);
        end
    end
    % Set to Extra
    Extra.qcov = C; Extra.muX = zeros(1,MCMCPar.n); Extra.invC = inv(C);
    % ---------------------------------------------------------------------

    Extra.save_in_memory = 'No';
    % Define data structures for use in computation of posterior density
    Measurement.MeasData = []; Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);
    % Define the example specific properties used to compute output
    Extra.mu = zeros(1,MCMCPar.n);          % Center point
    % What type of initial sampling
    Extra.InitPopulation = 'LHS_BASED';
    % Define modelName
    ModelName = 'normalfunc';
    % Define likelihood function
    option = 4;

end;

if example == 2,    % multivariate student t distribution with 60 degrees of freedom

    % Problem specific parameter settings
    MCMCPar.n = 25;                         % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 100000;                 % Maximum number of function evaluations
    MCMCPar.parallelUpdate = 0.9;           % Fraction of parallel direction updates

    % Recommended parameter settings
    MCMCPar.seq = 3;                        % Number of Markov Chains / sequences
    MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.nCR = 3;                        % Number of crossover values used
    MCMCPar.m0 = 10 * MCMCPar.n;            % Initial size of Z
    MCMCPar.k = 10;                         % Thinning parameter for appending X to Z
    MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
    MCMCPar.steps = 500;                    % Number of steps before calculating convergence diagnostics
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % --------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of probability of CR values
    % --------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage -------------------
    Extra.reduced_sample_collection = 'Yes'; %  Only collect thinned samples
    Extra.T = 1e1;                           % Every Tth sample is collected
    % --------------------------------------------------------------------------------------------

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [-5 * ones(1,MCMCPar.n)]; ParRange.maxn = [15 * ones(1,MCMCPar.n)];
    % Define the boundary handling
    Extra.BoundHandling = 'none';

    % ---------------------- Define covariance matrix ---------------------
    % Construct the dxd correlation matrix
    Extra.qcorr = 0.5*eye(MCMCPar.n) + 0.5*ones(MCMCPar.n);
    % ---------------------------------------------------------------------

    % How many degrees of freedom of student distribution used as target function?
    Extra.df = 60;

    % Make sure C is a valid covariance matrix
    [Extra.R,err] = cholcov(Extra.qcorr,0);

    % Define dimensionality
    Extra.d = MCMCPar.n;

    Extra.save_in_memory = 'No';
    % Define data structures for use in computation of posterior density
    Measurement.MeasData = []; Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);
    % Define the example specific properties used to compute output
    Extra.mu = zeros(1,MCMCPar.n);          % Center point
    % What type of initial sampling
    Extra.InitPopulation = 'LHS_BASED';
    % Define modelName
    ModelName = 'multi_student';
    % Define likelihood function
    option = 1;

end;

if example == 3, % n-dimensional banana shaped Gaussian distribution

    % Application specific parameter settings
    MCMCPar.n = 10;                         % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 50000;                  % Maximum number of function evaluations
    MCMCPar.parallelUpdate = 0.9;           % Fraction of parallel direction updates

    % Recommended parameter settings
    MCMCPar.seq = 3;                        % Number of Markov Chains / sequences
    MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.nCR = 3;                        % Number of crossover values used
    MCMCPar.m0 = 10 * MCMCPar.n;            % Initial size of Z
    MCMCPar.k = 10;                         % Thinning parameter for appending X to Z
    MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
    MCMCPar.steps = 100;                    % Number of steps before calculating convergence diagnostics
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % --------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of probability of CR values
    % --------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage -------------------
    Extra.reduced_sample_collection = 'Yes';% Only collect thinned samples
    Extra.T = 50;                           % Every Tth sample is collected
    % --------------------------------------------------------------------------------------------

    % Define the specific properties of the banana function
    Extra.mu   = [zeros(1,MCMCPar.n)];                      % Center of the banana function
    Extra.cmat = eye(MCMCPar.n); Extra.cmat(1,1) = 100;
    Extra.imat = inv(Extra.cmat);                           % Inverse of target covariance
    Extra.bpar = [0.1];                                     % "bananity" of the target, see bananafun.m

    % What type of initial sampling
    Extra.InitPopulation = 'COV_BASED';
    % Provide information to do alternative sampling
    Extra.muX = Extra.mu;                                   % Provide mean of initial sample
    Extra.qcov = eye(MCMCPar.n) * 5;                        % Initial covariance
    Extra.save_in_memory = 'No';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [-Inf * ones(1,MCMCPar.n)]; ParRange.maxn = [Inf * ones(1,MCMCPar.n)];

    % Define the boundary handling
    Extra.BoundHandling = 'None';
    % Define data structures for use in computation of posterior density
    Measurement.MeasData = []; Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);
    % Define modelName
    ModelName = 'Banshp';
    % Define likelihood function
    option = 4;

end;

if example == 4,    % n-dimensional multimodal mixture distribution

    % Problem specific parameter settings
    MCMCPar.n = 10;                         % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 200000;                 % Maximum number of function evaluations
    MCMCPar.parallelUpdate = 0.9;           % Fraction of parallel direction updates

    % Recommended parameter settings
    MCMCPar.seq = 5;                        % Number of Markov Chains / sequences
    MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.nCR = 3;                        % Number of crossover values used
    MCMCPar.m0 = 10 * MCMCPar.n;            % Initial size of Z
    MCMCPar.k = 10;                         % Thinning parameter for appending X to Z
    MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
    MCMCPar.steps = 100;                    % Number of steps before calculating convergence diagnostics
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % --------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of probability of CR values
    % --------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage -------------------
    Extra.reduced_sample_collection = 'Yes'; % Only collect thinned samples
    Extra.T = 10;                          % Every Tth sample is collected
    % --------------------------------------------------------------------------------------------

    % What type of initial sampling
    Extra.InitPopulation = 'COV_BASED';
    % Provide information to do alternative sampling
    Extra.muX = zeros(1,MCMCPar.n);         % Provide mean of initial sample
    Extra.qcov = eye(MCMCPar.n);            % Initial covariance

    Extra.Lam = eye(MCMCPar.n);             % covariance
    Extra.mu1 = -5 * ones(1,MCMCPar.n);     % center point of first density
    Extra.mu2 =  5 * ones(1,MCMCPar.n);     % center point of second density
    Extra.sigma = eye(MCMCPar.n);

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [-Inf * ones(1,MCMCPar.n)]; ParRange.maxn = [Inf * ones(1,MCMCPar.n)];
    % Define the boundary handling
    Extra.BoundHandling = 'None';
    % Save in memory or not
    Extra.save_in_memory = 'No';

    % Define data structures for use in computation of posterior density
    Measurement.MeasData = []; Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);
    % Define modelName
    ModelName = 'mixturemodel';
    % Define likelihood function
    option = 1;

end;

if example == 5,    % HYMOD rainfall - runoff model (coded in MATLAB)

    % Problem specific parameter settings
    MCMCPar.n = 5;                          % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 15000;                  % Maximum number of function evaluations
    MCMCPar.parallelUpdate = 0.9;           % Fraction of parallel direction updates
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % Recommended parameter settings
    MCMCPar.seq = 3;                        % Number of Markov Chains / sequences
    MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.nCR = 3;                        % Number of crossover values used
    MCMCPar.m0 = 10 * MCMCPar.n;            % Initial size of Z
    MCMCPar.k = 10;                         % Thinning parameter for appending X to Z
    MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
    MCMCPar.steps = 10;                     % Number of steps before calculating convergence diagnostics
 
    % --------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % --------------------------------------------------------------------------------------------
    
    % --------------------------------------- Added for reduced sample storage -------------------
    Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
    Extra.T = 1000;                         % Every Tth sample is collected
    % --------------------------------------------------------------------------------------------

    % What type of initial sampling
    Extra.InitPopulation = 'LHS_BASED';
    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [1.0 0.10 0.10 0.00 0.10]; ParRange.maxn = [500 2.00 0.99 0.10 0.99];
    % Define the boundary handling
    Extra.BoundHandling = 'Reflect';
    % Save in memory or not
    Extra.save_in_memory = 'Yes';

    % Load the Leaf River data
    load bound.txt;

    % Then read the boundary conditions -- only do two years
    Extra.MaxT = 795; 

    % Define the PET, Measured Streamflow and Precipitation.
    Extra.PET = bound(1:Extra.MaxT,5); Extra.Precip = sum(bound(1:Extra.MaxT,6:9),2);

    % Define the measured streamflow data
    Measurement.MeasData = bound(65:Extra.MaxT,4); Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);
    % Define modelName
    ModelName = 'hymodMATLAB';
    % Define likelihood function -- Sum of Squared Error
    option = 3;

end;

if example == 6,    % HYMOD rainfall - runoff model (coded in FORTRAN)

    % Problem specific parameter settings
    MCMCPar.n = 5;                          % Dimension of the problem (number of parameters to be estimated)
    MCMCPar.ndraw = 15000;                  % Maximum number of function evaluations
    MCMCPar.parallelUpdate = 0.9;           % Fraction of parallel direction updates
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % Recommended parameter settings
    MCMCPar.seq = 3;                        % Number of Markov Chains / sequences
    MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.nCR = 3;                        % Number of crossover values used
    MCMCPar.m0 = 10 * MCMCPar.n;            % Initial size of Z
    MCMCPar.k = 10;                         % Thinning parameter for appending X to Z
    MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
    MCMCPar.steps = 10;                     % Number of steps before calculating convergence diagnostics
 
    % --------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % --------------------------------------------------------------------------------------------
    
    % --------------------------------------- Added for reduced sample storage -------------------
    Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
    Extra.T = 1000;                         % Every Tth sample is collected
    % --------------------------------------------------------------------------------------------

    % What type of initial sampling
    Extra.InitPopulation = 'LHS_BASED';
    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [1.0 0.10 0.10 0.00 0.10]; ParRange.maxn = [500 2.00 0.99 0.10 0.99];
    % Define the boundary handling
    Extra.BoundHandling = 'Reflect';
    % Save in memory or not
    Extra.save_in_memory = 'Yes';

    % Leaf River data -- forcing conditions not needed --> externally loaded by FORTRAN executable
    load bound.txt; 
    
    % Then read the boundary conditions -- only do two years
    Extra.MaxT = 795; 

    % Define the measured streamflow data
    Measurement.MeasData = bound(65:Extra.MaxT,4); Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);

    % Define modelName
    ModelName = 'hymodFORTRAN';
    
    % Store working directory and subdirectory containing the files needed to run this example
	Extra.workdir = pwd; Extra.subdir = [pwd '\example_' num2str(example)];

    % Define likelihood function -- Sum of Squared Error
    option = 3;
end;

if example == 7,    % Rainfall-runoff model with generalized log-likelihood

    % ---------------------------- Check the following 2 papers ------------------------------- %
    %                                                                                           %
    % G. Schoups, J.A. Vrugt, F. Fenicia, and N.C. van de Giesen (2010), Corruption of          %
	%     accuracy and efficiency of Markov Chain Monte Carlo simulation by inaccurate          %
	%     numerical implementation of conceptual hydrologic models, Water Resources             %
    %     Research, 46, W10530, doi:10.1029/2009WR008648.                                       %
    %                                                                                           %
    % G. Schoups, and J.A. Vrugt (2010), A formal likelihood function for parameter and         %
    %     predictive inference of hydrologic models with correlated, heteroscedastic and        %
    %     non-Gaussian errors, Water Resources Research, 46, W10531, doi:10.1029/2009WR008933.  % 
    %                                                                                           %
    % ----------------------------------------------------------------------------------------- %

    % Problem specific parameter settings
    MCMCPar.n = 11;                         % Dimension of the problem
    MCMCPar.ndraw = 5000;                   % Maximum number of function evaluations
    MCMCPar.parallelUpdate = 0.9;           % Fraction of parallel direction updates
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % Recommended parameter settings
    MCMCPar.seq = 3;                        % Number of Markov Chains / sequences
    MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.nCR = 3;                        % Number of crossover values used
    MCMCPar.m0 = 10 * MCMCPar.n;            % Initial size of Z
    MCMCPar.k = 10;                         % Thinning parameter for appending X to Z
    MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
    MCMCPar.steps = 10;                     % Number of steps before calculating convergence diagnostics

    % --------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % --------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage -------------------
    Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
    Extra.T = 1000;                         % Every Tth sample is collected
    % --------------------------------------------------------------------------------------------

    % What type of initial sampling
    Extra.InitPopulation = 'LHS_BASED';
    % Give the parameter ranges (minimum and maximum values)
    %parno:       1     2     3     4     5     6     7      8    9     10    11   12   13   14   15   16   17   18   19   20  21
    %parname:     fA    Imax  Smax  Qsmax alE   alS   alF    Pf   Kfast Kslow std0 std1 beta xi   mu1  phi1 phi2 phi3 phi4 K   lambda
    Extra.fpar = [1     0     100   10    100   1e-6  1e-6   0    2     70    0.1  0    0    1    0    0    0    0    0    0   1];
    parmin =     [1     0     10    0     1e-6  1e-6 -10     0    0     0     0    0   -1    0.1  0    0    0    0    0    0   0.1 ];
    parmax =     [1     10    1000  100   100   1e-6  10     0    10    150   1    1    1    10   100  1    1    1    1    1   1];
    Extra.idx_vpar = [2 3 4 5 7 9 10 11 12 13 16];    
    ParRange.minn = parmin(Extra.idx_vpar); ParRange.maxn = parmax(Extra.idx_vpar);
    % Define the boundary handling
    Extra.BoundHandling = 'Reflect';
    % Save in memory or not
    Extra.save_in_memory = 'Yes';

    % Load the French Broad data
    daily_data = load('03451500.dly');

    % First two years are warm-up
    Extra.idx = [731:size(daily_data,1)]';

    % Define the PET, Measured Streamflow and Precipitation.
    Extra.Precip = daily_data(:,4);
    Extra.Ep     = daily_data(:,5);

    % Define the measured streamflow data
    Measurement.MeasData = daily_data(Extra.idx,6); 
    Measurement.Sigma = []; 
    Measurement.N = size(Measurement.MeasData,1);
    % Define modelName
    ModelName = 'hmodel';
    % Use generalized likelihood function
    option = 8; 
end;

if example == 8,	% HYDRUS-1D soil hydraulic model: using prior information on soil hydraulic parameters
	
	% -------------------------------- Check the following paper ------------------------------ %
    %                                                                                           %
    %   B. Scharnagl, J.A. Vrugt, H. Vereecken, and M. Herbst (2011), Bayesian inverse          % 
	%	modeling of soil water dynamics at the field scale: using prior information             % 
	%	on soil hydraulic properties, Hydrology and Earth System Sciences.                      %  
	%                                                                                           %
    % ----------------------------------------------------------------------------------------- % 

    % Problem specific parameter settings
    MCMCPar.n = 7;                          % Dimension of the problem
    MCMCPar.ndraw = 5000;                   % Maximum number of function evaluations
    MCMCPar.parallelUpdate = 0.9;           % Fraction of parallel direction updates
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % Recommended parameter settings
    MCMCPar.seq = 3;                        % Number of Markov Chains / sequences
    MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.nCR = 3;                        % Number of crossover values used
    MCMCPar.m0 = 10 * MCMCPar.n;            % Initial size of Z
    MCMCPar.k = 10;                         % Thinning parameter for appending X to Z
    MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
    MCMCPar.steps = 10;                     % Number of steps before calculating convergence diagnostics

    % --------------------------------------------------------------------------------------------
	Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % --------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage -------------------
    Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
    Extra.T = 1000;                         % Every Tth sample is collected
    % --------------------------------------------------------------------------------------------
	
	% Define the boundary handling
    Extra.BoundHandling = 'Reflect';
	
    % Save in memory or not
    Extra.save_in_memory = 'Yes';
	
	% Define feasible parameter space (minimum and maximum values)
	%				1		2		3				4			5			6		7     
	%				[thetar	thetas	log10(alpha)	log10(n)	log10(Ks)	L		hLB
	ParRange.minn =	[0.0430	0.4090	-2.5528			0.1790		-2.2366		-5.4900	-250];
	ParRange.maxn =	[0.0910 0.4810	-2.0706			0.2670		-0.0800		6.2700	-50];
	
	% Store working directory and subdirectory containing the files needed to run this example
	Extra.workdir = pwd; Extra.subdir = [pwd '\example_' num2str(example)];
	
	% Add subdirectory to search path
	addpath(Extra.subdir)
	
	% Provide observational data and data needed to modify the initial and boundary conditions
	[Measurement,Extra] = ProvideData(Extra);
	
	% Define the boundary handling
    Extra.BoundHandling = 'Reflect';

	% Define model name
	ModelName = 'HYDRUS';
	
	% Define option (model computes log-likelihood)
	option = 4;

	% Indicate the use prior information
	Extra.InitPopulation = 'PRIOR';
	
	% Specify the prior distributions for the various parameters
	Extra.prior = {'normrnd(0.0670,0.0060)',...
				   'normrnd(0.4450,0.0090)',...
				   'normrnd(-2.310,0.0600)',...
				   'normrnd(0.2230,0.0110)',...
				   'normrnd(-1.160,0.2700)',...
				   'normrnd(0.3900,1.4700)',...
				   'unifrnd(-250,-50)'};
	
end;

% No restart -- just running for the first time. Restart can be used if run
% is termined while running; Use Restart = 'Yes'; to reinitialize the code
% from the saved files. 
Restart = 'No';

% Run the distributed DREAM algorithm with sampling from past
[Sequences,Reduced_Seq,X,Z,output] = dream_zs(MCMCPar,ParRange,Measurement,ModelName,Extra,option,Restart);