clear all
close all
clc
%% initialisation des variables
global sigma epsilon k0 Natome m

m = 1 ;	% masse d'un atome en kg
sigma = 1;  % distance ou le potentiel s'annule en m
epsilon = 1;    % porfonderur du puit du potentiel
dt = 0.005;	% pas du temps en s
Niter = 3e4;	% nombre d'iterations
kB = 1;	% cts de boltzman
T = 3;	% temperature en K
k0=30;  % raideur 
ftrac=[0 0 20];	% force de traction en N
% df=15/Niter;
Natome=20;   % nbr d'atome au milieu

%% initialisation des positions
P=zeros(Natome+1,3,Niter+1);    % positions des particules
P(:,:,1)=[zeros(Natome+1,1) zeros(Natome+1,1) 1.5*(1:Natome+1)'];
% P(:,:,1)=[0 0 0; 1 2 3];
% sigma1 = zeros (Niter,1) ; 
%% calcul de vitesse aleatoire
v=(3*kB*T/m)^0.5;
Vi=vitesse(v,Natome+1);
Vi=cancelTrans(Vi);
Vi=cancelRot(P(:,:,1),Vi);
%% rescaling des vitesses
vt=mean(sqrt(sum(Vi.^2,2)),1);
Tt=m*vt^2/3/kB;
lambda=(Tt/T)^0.5;
Vi=Vi/lambda;
vt=mean(sqrt(sum(Vi.^2,2)),1);
V=zeros(Natome+1,3,Niter+1);
V(:,:,1)=Vi;

%% video
writerObj=VideoWriter('demorand111');
open(writerObj);
figure('numbertitle','off','name','iteration des atomes dans un thermostat')

%%
timeflag=1e4;   % nombre d'iterations d'evolution libre
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
%     
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
%     
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

    sigma1(i,1)=sqrt(sum(forcetot(P(Natome,:,i),P(Natome+1,:,i),P(:,:,i),Natome+1).^2,2));
    
%     video
	if(mod(i,25)==0)
    	Ptemp=[zeros(1,3);P(:,:,i)]; 
%         zeros(1,3);
            if(i<timeflag)
                plot3(Ptemp(:,1),Ptemp(:,2),Ptemp(:,3),'.-r','MarkerSize',25);
            end
            if(i>=timeflag)
                plot3(Ptemp(:,1),Ptemp(:,2),Ptemp(:,3),'.-b','MarkerSize',25);
            end
        title('evolution des atomes dans un thermostat');xlabel('x(sigma)');ylabel('y(sigma)');zlabel('z(sigma)')
%         view(0,0);
        axis([-20 20 -20 20 -20 60]);
        grid
        
        frame=getframe(gcf);
        frame.cdata = imresize(frame.cdata, [500 800]);
        writeVideo(writerObj,frame);
        drawnow;
	end
end
close(writerObj);
% for i=1:Niter
%     Def(i)=(sqrt(sum(P(Natome+1,:,i).^2,2))-sqrt(sum(P(Natome+1,:,1).^2,2)))/sqrt(sum(P(Natome+1,:,1).^2,2));
% %     Def(i)=sqrt(sum(P(Natome+1,:,i).^2,2)
% %     Def(i)=sqrt(sum(P(Natome+1,:,i).^2,2))-sqrt(sum(P(Natome+1,:,1).^2,2));
% end
% sigma1(1,1)=sqrt(sum(forcetot(P(Natome,:,1),P(Natome+1,:,1),P(:,:,1),Natome+1).^2,2));
% figure(2)
% C=polyfit(Def,5+dt*[1:size(Def,2)],1);
% 
% fplot(@(x) C(1)*x+C(2),[0 2])
% plot(Def',sigma1,'r')
% defmoy=mean(Def)
% L1=polyfit(Def,sigma1,1) ;
% Eyoung1=L1(1)

% for i=1:Niter+1
% vtt(i)=mean(sqrt(sum(V(1:Natome+1,:,i).^2,2)),1);
% Ttt(i)=m*vtt(i)^2/3/kB;
% end
% figure('numbertitle','off','name','variation de la temperature dans un thermostat')
% plot(Ttt);
% title('variation de la temperature dans un thermostat');xlabel('iteration');ylabel('temperature (K)');
% legend('temperature moyenne des atomes','Location','northeast')
% saveas(gcf,'Temperature elastomere.jpg')

% figure(4)
% plot(vtt);




