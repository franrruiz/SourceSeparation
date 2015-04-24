%function ACC_REDD(H,Q)
Nd=5;
T=1354;
Tsubsample=250;
BASEDIR1=['PCCdata16kHz_isolated/resultsFFBS_lap/S' num2str(Nd) '_T' num2str(T) '_Tsub' num2str(Tsubsample)];
load([BASEDIR1 '/Final.mat'],'data','init','samples','samplesAll', 'LLH', 'M_EST');

ACC=zeros(1,2000);
cad_ord=cell(1,2000);
devices=double(data.speakers~=0)';
for it=1:2000
    cadenas=double(samplesAll{it}.Z~=0);
    [ACC(it) cad_ord{it}]= calculaAccuracy(cadenas,devices);
end

BASEDIR1=['PCCdata16kHz_isolated/resultsFFBS_lap'];
save([BASEDIR1 '/ADER_S' num2str(Nd) '_T' num2str(T) '_Tsub' num2str(Tsubsample) '.mat'],'ACC','cad_ord','devices');