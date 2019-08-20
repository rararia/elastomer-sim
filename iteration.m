function defmoy = iteration(T,ftrac,videoflag,Pini)
global sigma epsilon k0 Natome m dt Niter kB
%% initialisation des positions
P=zeros(Natome+1,3,Niter+1);    % positions des particules
P(:,:,1)=Pini;  % positions initiales
PosG=zeros(Niter+1,3);  % position du centre de masse
PosG(1,:)=centerDeMasse(Pini);
sigma1 = zeros (Niter,1) ;  % force de rappel
%% calcul de vitesse aleatoire
v=(3*kB*T/m)^0.5;   % norme de vitesse
Vi=vitesse(v,Natome+1); % initialisation des atomes
Vi=cancelTrans(Vi);
Vi=cancelRot(P(:,:,1),Vi);
%% rescaling des vitesses
vt=mean(sqrt(sum(Vi.^2,2)),1);
Tt=m*vt^2/3/kB;
lambda=(Tt/T)^0.5;  % coeff de changement
Vi=Vi/lambda;   % rescaling
vt=mean(sqrt(sum(Vi.^2,2)),1);
V=zeros(Natome+1,3,Niter+1);    % vitesses des particules
V(:,:,1)=Vi;    % initialisation
%%
timeflag=3e4;   % nombre d'iterations d'evolution libre
df=ftrac/(Niter-timeflag);  % pas de la force
ftemp=[0 0 0];
figure(1)
%% iteration pour tous
% schema de Verlet avec un pas
% P(n+1)=P(n)+dt*V(n)+dt^2/2*F(n)
% V(n+1)=V(n)+dt/2*(F(n+1)+F(n))
for i=1:Niter
    if(i>=timeflag)
        ftemp=ftemp+df; % excercer la force de traction
    end
%% calcul de P
    for k=2:Natome
        fextcurrent=forcetot(P(k+1,:,i),P(k,:,i),P(:,:,i),k)+forcetot(P(k-1,:,i),P(k,:,i),P(:,:,i),k);
        P(k,:,i+1)=P(k,:,i)+dt*V(k,:,i)+dt^2/2*fextcurrent/m;
%         V(k,:,i+1)=V(k,:,i)+dt*fextcurrent/m;        
    end
    fext1current=forcetot([0,0,0],P(1,:,i),P(:,:,i),1)+forcetot(P(2,:,i),P(1,:,i),P(:,:,i),1);
    P(1,:,i+1)=P(1,:,i)+dt*V(1,:,i)+dt^2/2*fext1current/m;

%     V(1,:,i+1)=V(1,:,i)+dt/2*fext1current/m;
    
    fextdercurrent=forcetot(P(Natome,:,i),P(Natome+1,:,i),P(:,:,i),Natome+1)+ftemp;
    P(Natome+1,:,i+1)=P(Natome+1,:,i)+dt*V(Natome+1,:,i)+dt^2/2*fextdercurrent/m;
%     V(Natome+1,:,i+1)=V(Natome+1,:,i)+dt*fextdercurrent/m;

%% calcul de V
    for k=2:Natome
        fextcurrent=forcetot(P(k+1,:,i),P(k,:,i),P(:,:,i),k)+forcetot(P(k-1,:,i),P(k,:,i),P(:,:,i),k);
        fextnext=forcetot(P(k+1,:,i+1),P(k,:,i+1),P(:,:,i+1),k)+forcetot(P(k-1,:,i+1),P(k,:,i+1),P(:,:,i+1),k);
        V(k,:,i+1)=V(k,:,i)+dt/2*(fextcurrent/m+fextnext/m);        
    end
    fext1next=forcetot([0,0,0],P(1,:,i+1),P(:,:,i+1),1)+forcetot(P(2,:,i+1),P(1,:,i+1),P(:,:,i+1),1);
    V(1,:,i+1)=V(1,:,i)+dt/2*(fext1current/m+fext1next/m);    
    
    fextdernext=forcetot(P(Natome,:,i+1),P(Natome+1,:,i+1),P(:,:,i+1),Natome+1)+ftemp;
    V(Natome+1,:,i+1)=V(Natome+1,:,i)+dt/2*(fextdercurrent/m+fextdernext/m);
    
%     ftrac(3)=ftrac(3)+df;
%%    
%     thermosta
    vi=mean(sqrt(sum(V(1:Natome+1,:,i+1).^2,2)),1);
    Ti=m*vi^2/3/kB;
    lambda=sqrt(Ti/T);
    V(1:Natome+1,:,i+1)=V(1:Natome+1,:,i+1)/lambda;
%     force de rappel
%     sigma1(i,1)=sqrt(sum(forcetot(P(Natome,:,i),P(Natome+1,:,i),P(:,:,i),Natome+1).^2,2));    
    
%     video
    PosG(i+1,:)=centerDeMasse([zeros(1,3);P(:,:,i)]);
    if(videoflag~=0)
        if(mod(i,100)==0)
            Ptemp=[zeros(1,3);P(:,:,i)];  
%             subplot(1,2,1)
%             change color
            if(i<timeflag)
                plot3(Ptemp(:,1),Ptemp(:,2),Ptemp(:,3),'.-r','MarkerSize',25);
            end
            if(i>=timeflag)
                plot3(Ptemp(:,1),Ptemp(:,2),Ptemp(:,3),'.-b','MarkerSize',25);
            end
            axis([-10 10 -10 10 -40 60]);
            grid
%            view(0,0);
%             evolution du centre de masse
%             subplot(1,2,2)
%             plot3(PosG(i+1,1),PosG(i+1,2),PosG(i+1,3),'.r','MarkerSize',25);
%             axis([-10 10 -10 10 -40 60]);
%             grid
            drawnow;
        end
    end
    
end
Def=zeros((Niter-timeflag),1);
balance=0;
for i=(timeflag-1.5e4):timeflag
    balance=balance+P(Natome+1,3,i+timeflag);   % point balance
end
balance=balance/(timeflag-1.5e4);

% balance=P(Natome+1,3,timeflag);

% balance=0;
if(sum(ftrac.^2)~=0)
    for i=1:(Niter-timeflag)
        Def(i)=P(Natome+1,3,i+timeflag)-balance;
%         Def(i)=(sqrt(sum(P(Natome+1,:,i).^2,2))-sqrt(sum(P(Natome+1,:,1).^2,2)))/sqrt(sum(P(Natome+1,:,1).^2,2));
%         Def(i)=sqrt(sum(P(Natome+1,:,i).^2,2);
%         Def(i)=sqrt(sum(P(Natome+1,:,i).^2,2))-sqrt(sum(P(Natome+1,:,1).^2,2));
    end
end
fff=sqrt(sum(ftrac.^2,2));  % norme de force max
ffft=[1:(Niter-timeflag)]*fff/(Niter-timeflag); % evolution de force
figure(1)
plot(Def,ffft)
title(['deformation de ' num2str(Natome) '  atomes dans un thermostat avec T=' num2str(T) 'K']);xlabel('def z(sigma)');ylabel('force(epsilon/sigma)');
saveas(gcf,[num2str(Natome) 'atomes, T = ' num2str(T) '.jpg'])
defmoy=mean(Def);
% posgmoy=mean(PosG,1)
% if(sum(ftrac.^2)~=0)    
%     L1=polyfit(Def,sigma1,1) ;
%     Eyoung=L1(1)
% end
end

