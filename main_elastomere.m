clear all
close all
clc
%% initialisation des variables
global sigma epsilon k0 Natome m dt Niter kB

m = 1 ;	% masse d'un atome en kg
sigma = 1;  % distance ou le potentiel s'annule en m
epsilon = 1;    % profonderur du puit du potentiel
dt = 0.002;	% pas du temps en s
Niter = 10e4;	% nombre d'iterations
kB = 1;	% cts de boltzman
T = 5;	% temperature en K
k0=30;  % raideur 
ftrac=[0 0 0];	% force de traction en N

Natome=40;   % nbr d'atome au milieu

videoflag=1;    % indicateur si on fait la video ou pas
Pini=[zeros(Natome+1,1) zeros(Natome+1,1) 1.5*(1:Natome+1)'];   % position initiale

%% demo libre
% defmoy=iteration(T,ftrac,videoflag,Pini);

%% traction avec une force qui augmente lineairement
Def=zeros(5,1);
videoflag=0;
% figure(1)
tempt=[1,2,3,5,10]; % temperatures
for i=1:size(tempt,2)
    T=tempt(i);
%     for j=1:size(Def,2)
        ftrac=[0 0 20]; % force de traction maximum
        Def(i)=iteration(T,ftrac,videoflag,Pini);     
%     end
%     plot(Def(i,:),20*(1:size(Def,2)));
%     title('deformation des atomes dans un thermostat');xlabel('def relative');ylabel('force(N)');
% 
%     axis([0 30 0 5]);
%     hold on
end
%     legend('T=2.5','T=5','T=7.5','Location','northeast')
%     saveas(gcf,'Temperature elastomere.jpg')