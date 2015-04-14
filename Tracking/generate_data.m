function generate_data(Nt,s2y,T)
Ts=0.01;
Ns=25;
pathL=1.8;
W=20;
s2u=1;
d0=1;
Pt=1;
Gx=[1 0 Ts 0; 0 1 0 Ts; 0 0 1 0; 0 0 0 1];
Gu= [Ts^2/2 0; 0 Ts^2/2; Ts 0; 0 Ts];

pos=0:5:W;
[a b] = meshgrid(pos);
sensors=[a(: ) b(:)];

state=zeros(Nt,4,T);
for nt=1:Nt
    Tini = randi([1 round(T/2)],1,1);
    Tend = min(Tini+round(T/2)-1,T);
    state (nt,:,Tini)= [W*rand(2,1); randn(2,1)];
    for t=Tini+1:Tend
       state (nt,:,t)=Gx*squeeze(state (nt,:,t-1))'+sqrt(s2u)*Gu*randn(2,1);
    end
end

Ptx=zeros(Ns,T);
for t=1:T
    d=zeros(Ns,1);
    for nt=1:Nt
       if state (nt,1,t)~=0
            d= d+1./(((sensors(:,1)-state (nt,1,t)).^2 +(sensors(:,2)-state (nt,2,t)).^2).^(pathL/2));
       end
    end
    if d>0
        Ptx(:,t)= 10^(Pt/10)*d0^pathL*d;
    end
end

obs=Ptx+sqrt(s2y)*randn(Ns,T);

save (['dataTracking_Nt' num2str(Nt) '_s2y' num2str(s2y) '_T' num2str(T) '.mat']);
