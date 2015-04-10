function sim_Tracking(noiseVar,T,Nd,Ns,Niter,LastIt)

addpath(genpath('sampleFunc/'));
addpath(genpath('auxFunc/'));

randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

flagRecovered=0;
if LastIt>0
    flagRecovered=1;
end

%% Configuration parameters
param.Nd = Nd;                        % Number of devices
param.D = Ns;                        %Dimensionality of the observations = number of sensors
param.T  = T;                         % Length of the sequence
param.L = 1;
param.flag0 = 1;    % Consider symbol 0 as part of the constellation (if false, transmitters are always active)
param.Niter = Niter;  % Number of iterations of the sampler
param.saveCycle = 200;
param.storeIters = 2000;
param.constellation=1;
param.d0=1;
param.pathL= 1.8;
param.Ts=0.01;

%% Generate data
BASEDIR1=['Tracking/resultsPGAS/M' num2str(param.Nd) '_T' num2str(T) '_s2y' num2str(noiseVar)];
if(~isdir(BASEDIR1))
    mkdir(BASEDIR1);
end
data.s2y=noiseVar;
data.s2u=1;
data.Ptx= 0.01; % Transmitted power in dB
data.Gx=[1 0 param.Ts 0; 0 1 0 param.Ts; 0 0 1 0; 0 0 0 1];
data.Gu= [param.Ts^2/2 0; 0 param.Ts^2/2; param.Ts 0; 0 param.Ts];
data.W=200;

[data.obs  data.states data.sensors]= generate_data(param.T, param.Nd, param.D, param.d0, param.pathL,data.s2y,data.s2u,data.Ptx, data.W, data.Gx, data.Gu);

%% Configuration parameters for BCJR, PGAS, EP, FFBS and collapsed Gibbs
param.bcjr.p1 = 0.95;
param.bcjr.p2 = 0.05;
param.pgas.N_PF = 3000;
param.pgas.N_PG = 3000;
param.pgas.Niter = 1;
param.pgas.returnNsamples = 1;
param.pgas.maxM = 40;
%param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');
param.ep.eps = 5e-7;
param.ep.beta = 0.2;
param.ep.Niter = 15;
param.colGibbs.Niter = 1;
param.ffbs.Niter = 1;

%% Configuration parameters for BNP and inference method
param.infer.symbolMethod = 'pgas';
param.infer.sampleNoiseVar = 0;
param.bnp.betaSlice1 = 0.5;
param.bnp.betaSlice2 = 5;
param.bnp.maxMnew = 15;
param.bnp.Mini = 1;


%% Hyperparameters
hyper.alpha = 1;    % Concentration parameter for Z ~ IBP(alpha)
hyper.gamma1 = 0.1; % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.gamma2 = 2;   % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.tau = 1;      % Parameter for s2y ~ IG(tau,nu)
hyper.nu = 1;       % Parameter for s2y ~ IG(tau,nu)

%% Initialization
if(~flagRecovered)
    if param.infer.sampleNoiseVar
        init.s2y = 20*rand;      % INITIALIZE s2y TO THE GROUND TRUTH
    else
        init.s2y = data.s2y;
    end
    init.am = 0.95*ones(param.bnp.Mini,1);
    init.bm = 0.05*ones(param.bnp.Mini,1);
    init.Z = zeros(param.bnp.Mini,4,param.T);
    init.nest = zeros(2,2,param.bnp.Mini);
    init.nest(1,1,:) = param.T;
    init.slice = 0;
    init.epAcc = 0;
    samples = init;
    samplesAll = cell(1,param.storeIters);
else
    load([BASEDIR1 '/it' num2str(LastIt) '.mat'],'data','init','samples','samplesAll');
    
end

%% Inference
for it=LastIt+1:param.Niter
    %% Algorithm
    
    % Step 1)
    % -Sample the slice variable
    samples.slice = sample_post_slice(data,samples,hyper,param);
    % -Sample new sticks (and the corresponding new parameters)
    samples = sample_newsticks(data,samples,hyper,param);
    
    % Step 2)
    % -Sample the symbols Z
    [samples.Z samples.seq samples.nest out] = sample_post_Z(data,samples,hyper,param);
    % -Compute some statistics of interest
    if(strcmp(param.infer.symbolMethod,'pgas'))
        
    elseif(strcmp(param.infer.symbolMethod,'ep'))
        samples.epAcc = samples.epAcc+out;
    end
    
    % Step 3)
    % -Remove unused chains
    samples = sample_remove_unused(data,samples,hyper,param);
    
    % Step 4)
    % -Sample the transition probabilities (semi-ordered construction)
    [samples.am samples.bm]= sample_post_transitionProb(data,samples,hyper,param);
    
    % Step 5)
    % -Sample the noise variance
    samples.s2y = sample_post_s2y(data,samples,hyper,param);

    
    %% Evaluation
%     % Trace of the estimated number of transmitters
     M_EST(it) = sum(sum(samples.Z(:,1,:)~=0,3)>0);
%     % Trace of the log-likelihood
     LLH(it) = 0;%compute_llh(data,samples,hyper,param);
    
    %% Save temporary result file
    if(mod(it,param.saveCycle)==0)
        save([BASEDIR1 '/it' num2str(it) '.mat'],'data','init','samples','samplesAll', 'LLH', 'M_EST');
        % If successfully saved, delete previous temporary file
%         if(exist([saveFile '/it' num2str(it-param.saveCycle) '.mat'],'file'))
%             delete([saveFile '/it' num2str(it-param.saveCycle) '.mat']);
%         end
    end
end
save([BASEDIR1 '/Final.mat'],'data','init','samples','samplesAll', 'LLH', 'M_EST');


