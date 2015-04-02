function s2y = sample_post_s2y(data,samples,hyper,param)

if(~param.infer.sampleNoiseVar)
    s2y = samples.s2y;
    return;
end

% Obtain parameters from the structs
Nr = param.D;
[Nt aux T] = size(samples.Z);
L = param.L;

% Build matrix S containing the symbols and their shifted replicas,
% and matrix H containing the channel coefficients
Ptot=zeros(Nr,T);
for itm=1:Nt 
    d= sqrt((repmat(data.sensors(:,1),1,T)-repmat(squeeze(samples.Z(itm,1,:))',Nr,1)).^2 +(repmat(data.sensors(:,2),1,T)-repmat(squeeze(samples.Z(itm,2,:))',Nr,1)).^2);
    Ptot= Ptot-param.pathL*log10(d);
end    
Ptot=Ptot+data.Ptx+param.pathL*log(param.d0);
% Posterior parameters
nuP = hyper.nu+T*Nr;
tauP = hyper.tau+sum(sum(abs(data.obs-Ptot).^2));
% [Note: There is a typo in the ModelNotes.pdf file (in that file, there is
%  an extra 0.5 factor for tauP). Here it has been fixed]

% Sample the noise variance from an Inverse-Gamma distribution
s2y = 1/gamrnd(nuP,1/tauP);

