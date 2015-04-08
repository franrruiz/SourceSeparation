function [obs state sensors]=generate_data(T, Nt, Ns, d0, pathL, s2y,s2u,Pt, W, Gx, Gu)
sensors= W.*rand(Ns,2);

state=zeros(Nt,4,T);
Ptx=zeros(Ns,T);
for nt=1:Nt
    Tini = randi([1 round(T/2)],1,1);
    Tend = min(Tini+round(T/2)-1,T);
    state (nt,:,Tini)= [W*rand(2,1); randn(2,1)];
    for t=Tini+1:Tend
       state (nt,:,t)=Gx*squeeze(state (nt,:,t-1))'+sqrt(s2u)*Gu*randn(2,1);
    end
end

for t=1:T
    d=zeros(Ns,1);
    for nt=1:Nt
       if state (nt,1,t)~=0
            d= d+sqrt((sensors(:,1)-state (nt,1,t)).^2 +(sensors(:,2)-state (nt,2,t)).^2);
       end
    end
    if d>0
        Ptx(:,t)= Pt+10*pathL*log10(d0)-10*pathL*log10(d);
    end
end

obs=Ptx+sqrt(s2y)*randn(Ns,T);
