%% initialisation des variables
global sigma epsilon k0 Natome m

m = 1 ;             % masse d'un atome en kg
sigma = 1;          % distance ou le potentiel s'annule en m
epsilon = 1;        % porfonderur du puit du potentiel
dt = 0.005;         % pas du temps en s
Niter = 2e4;        % nombre d'iterations
kB = 1;             % cts de boltzman
T = 1;              % temperature en K
k0=30;              % raideur 
vtrac=[0 0 0.1];    % vitesse de traction en m/s
dx=vtrac*dt;        % deplacement par un pas du temps en m
Natome=20;          % nbr d'atome au milieu
video=0;            % video tag
%% initialisation des positions
P1=zeros(Natome+1,3,Niter+1);                                           % positions pour etat libre
P2=zeros(Natome+1,3,Niter+1);                                           % positions pour le tirage
P1(:,:,1)=[zeros(Natome+1,1) zeros(Natome+1,1) 1.5*(1:Natome+1)'];
sigma1 = zeros (Niter,1) ;                                              % force de rappel
%% calcul de vitesse aleatoire
v=(3*kB*T/m)^0.5;
Vi=vitesse(v,Natome);
% suppression de la translation et la rotation
Vi=cancelTrans(Vi);
Vi=cancelRot(P1(:,:,1),Vi);
%% rescaling des vitesses
vt=mean(sqrt(sum(Vi.^2,2)),1);      
Tt=m*vt^2/3/kB;
lambda=(Tt/T)^0.5;
Vi=Vi/lambda;
V1=zeros(Natome+1,3,Niter+1);
V2=zeros(Natome+1,3,Niter+1);
V1(:,:,1)=[Vi;vtrac];
%% iteration pour etat libre
for i=1:Niter
    %calculer les positions
    for k=2:Natome
        fextcurrent=forcetot(P1(k+1,:,i),P1(k,:,i),P1(:,:,i),k)+forcetot(P1(k-1,:,i),P1(k,:,i),P1(:,:,i),k);
        P1(k,:,i+1)=P1(k,:,i)+dt*V1(k,:,i)+dt^2/2*fextcurrent/m;
    end
    fext1current=forcetot([0,0,0],P1(1,:,i),P1(:,:,i),1)+forcetot(P1(2,:,i),P1(1,:,i),P1(:,:,i),1);
    P1(1,:,i+1)=P1(1,:,i)+dt*V1(1,:,i)+dt^2/2*fext1current/m;
    fextdercurrent=forcetot(P1(Natome,:,i),P1(Natome+1,:,i),P1(:,:,i),Natome+1);
    P1(Natome+1,:,i+1)=P1(Natome+1,:,i)+dt*V1(Natome+1,:,i)+dt^2/2*fextdercurrent/m;

%calculer les vitesses
     for k=2:Natome
        fextcurrent=forcetot(P1(k+1,:,i),P1(k,:,i),P1(:,:,i),k)+forcetot(P1(k-1,:,i),P1(k,:,i),P1(:,:,i),k);
        fextnext=forcetot(P1(k+1,:,i+1),P1(k,:,i+1),P1(:,:,i+1),k)+forcetot(P1(k-1,:,i+1),P1(k,:,i+1),P1(:,:,i+1),k);
        V1(k,:,i+1)=V1(k,:,i)+dt/2*(fextcurrent/m+fextnext/m);
     end
     fext1next=forcetot([0,0,0],P1(1,:,i+1),P1(:,:,i+1),1)+forcetot(P1(2,:,i+1),P1(1,:,i+1),P1(:,:,i+1),1);
     V1(1,:,i+1)=V1(1,:,i)+dt/2*(fext1current/m+fext1next/m);
     fextdernext=forcetot(P1(Natome,:,i+1),P1(Natome+1,:,i+1),P1(:,:,i+1),Natome+1);
     V1(Natome+1,:,i+1)=V1(Natome+1,:,i)+dt/2*(fextdercurrent/m+fextdernext/m);
    
%     thermostat
    vi=mean(sqrt(sum(V1(1:Natome,:,i+1).^2,2)),1);
    Ti=m*vi^2/3/kB;
    lambda=sqrt(Ti/T);
    V1(1:Natome,:,i+1)=V1(1:Natome,:,i+1)/lambda;
%     video
if(video)
    if(mod(i,100)==0)
    	Ptemp=[zeros(1,3);P1(:,:,i)];  
        plot3(Ptemp(:,1),Ptemp(:,2),Ptemp(:,3),'.-r','MarkerSize',25);
        title('etat libre dans un thermostat');
        axis([-5 5 -5 5 -10 100]);
        grid
        drawnow;
    end
end
end

P2(:,:,1)=P1(:,:,Niter+1);
V2(:,:,1)=V1(:,:,Niter+1);
for i=1:Niter+1
    P2(Natome+1,:,i)=P2(Natome+1,:,1)+dx*(i-1);
    V2(Natome+1,:,i)=vtrac; 
end
%%
%iteration pour le tirage
for i=1:Niter
    %calculer les positions
    for k=2:Natome
        fextcurrent=forcetot(P2(k+1,:,i),P2(k,:,i),P2(:,:,i),k)+forcetot(P2(k-1,:,i),P2(k,:,i),P2(:,:,i),k);
        P2(k,:,i+1)=P2(k,:,i)+dt*V2(k,:,i)+dt^2/2*fextcurrent/m;
    end
    fext1current=forcetot([0,0,0],P2(1,:,i),P2(:,:,i),1)+forcetot(P2(2,:,i),P2(1,:,i),P2(:,:,i),1);
    P2(1,:,i+1)=P2(1,:,i)+dt*V2(1,:,i)+dt^2/2*fext1current/m;
%calculer les vitessses
     for k=2:Natome
        fextcurrent=forcetot(P2(k+1,:,i),P2(k,:,i),P2(:,:,i),k)+forcetot(P2(k-1,:,i),P2(k,:,i),P2(:,:,i),k);
        fextnext=forcetot(P2(k+1,:,i+1),P2(k,:,i+1),P2(:,:,i+1),k)+forcetot(P2(k-1,:,i+1),P2(k,:,i+1),P2(:,:,i+1),k);
        V2(k,:,i+1)=V2(k,:,i)+dt/2*(fextcurrent/m+fextnext/m);
     end
     fext1next=forcetot([0,0,0],P2(1,:,i+1),P2(:,:,i+1),1)+forcetot(P2(2,:,i+1),P2(1,:,i+1),P2(:,:,i+1),1);
     V2(1,:,i+1)=V2(1,:,i)+dt/2*(fext1current/m+fext1next/m);
    
%     thermostat
    vi=mean(sqrt(sum(V2(1:Natome,:,i+1).^2,2)),1);
    Ti=m*vi^2/3/kB;
    lambda=sqrt(Ti/T);
    V2(1:Natome,:,i+1)=V2(1:Natome,:,i+1)/lambda;
    
    sigma1(i)=sqrt(sum(forcetot(P2(Natome,:,i),P2(Natome+1,:,i),P2(:,:,i),Natome+1).^2,2));

%     video
if(video)
    if(mod(i,100)==0)
    	Ptemp=[zeros(1,3);P2(:,:,i)];  
        plot3(Ptemp(:,1),Ptemp(:,2),Ptemp(:,3),'.-b','MarkerSize',25);
        title('le tirage delastomere dans un thermostat');
        axis([-5 5 -5 5 -10 100]);
        grid
        drawnow;
    end
end
end

%%
%calcul de la deformation
for i=1:Niter
    Def(i)=(sqrt(sum(P2(Natome+1,:,i).^2,2))-sqrt(sum(P2(Natome+1,:,1).^2,2)))/sqrt(sum(P2(Natome+1,:,1).^2,2));
end
sigma1(1,1)=sqrt(sum(forcetot(P2(Natome,:,1),P2(Natome+1,:,1),P2(:,:,1),Natome+1).^2,2));
figure(2)
plot(Def,sigma1);
grid on
title('la force de rappel en fonction de la deformation(sans moyenne glissante)')
xlabel('la deformation relative')
ylabel('la force de rappel')

%%
% moyenne glissante(intervalle de 100, glisser 10 chaque fois)
Defm=zeros(990,1);
sigmam=zeros(990,1);
for i=1:10:9900
Defm((i-1)/10+1)=mean(Def(i:i+100));
sigmam((i-1)/10+1)=mean(sigma1(i:i+100));
end
%%
figure(3)
plot(Defm,sigmam);
grid on
title('la force de rappel en fonction de la deformation(avec moyenne glissante)')
xlabel('la deformation relative')
ylabel('la force de rappel')
%calcul le pente(module de Young)
L1=polyfit(Defm,sigmam,1) ;
Eyoung1=L1(1)
%verification du thermostat
for i=1:Niter+1
vtt(i)=mean(sqrt(sum(V2(1:Natome,:,i).^2,2)),1);
Ttt(i)=m*vtt(i)^2/3/kB;
end
figure(4)
plot(Ttt);





