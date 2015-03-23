close all

speakers=zeros(60*16000,2,5);

for sp=1:5
    t=0;
    for sen=1:4
        [data f]=wavread(['data/s' num2str(sp) '_' num2str(sen) '.wav']);
        t=t+poissrnd((20.*rand)*16000);
        speakers(t:t+size(data,1)-1,:,sp)=data;
        t=t+size(data,1)-1;
    end
end

 plot(squeeze(speakers(1:500:end,1,:)))
 
 save('data2.mat','speakers')