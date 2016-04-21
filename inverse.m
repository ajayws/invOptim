clc;clear all;close all

% Random generator. Hapus bila tidak memiliki preferensi
% load('PSO-1.mat')
% rng(S)
% S=rng; % Generate new random seed. Gunakan ini bila tidak memiliki PSO-1.mat

% Parameter model sintetik
a=20;
xobs=[5:5:45]'; % reciever linear
% xobs=[25];
A=1:50;

% Parameter pendekatan linear
epsilon=2.5e-7; % Faktor redaman. Perbesar nilainya bila inversi tidak stabil
e=5e-1; % ambang batas error diterima.
maxiter=1000; %iterasi maksimum
deltaM=3; % dalam kilometer.
M0=[3]; % dugaan awal

% Parameter PSO
c1=0.3; 
c2=0.7;
w=0.1;
Nagent=15; % Nilai agent yang diinginkan
Apso=rand(1,Nagent)*50; % randomisasi posisi awal agent

% Parameter GA
maxiterga=100;
Nbin=6;
co=3;
Popset=5;
Npop=Popset*3+1;
Aga=rand(1,Npop)*50;

%Parameter Simulated annealing
T0=50;
T=T0;

%% forward modelling
yobs=linear1d(a,0,xobs);
noise=1*(rand(size(yobs))-0.5);
yobs=yobs+noise;

%% grid search
fprintf('Metode Grid Search :\n')
tic
for i=1:size(A,2)
    yi=linear1d(A(i),0,xobs);
    t(i)=sqrt(1/numel(yobs)*sum((yi-yobs).^2));
end
waktu=toc;
[Ma,I] = min(t(:));
% Tampilan di console
fprintf('Hasil :\n')
fprintf('    a : %2i\n',A(I))
fprintf('Waktu running : %7.5fs\n\n',waktu)
% Plot Spasial
figure
plot(A,t)
hold on
plot(a,0,'*r')
title('Grid Search')
ylabel('missfit'),xlabel('model parameter a')
axis([0 50 0 1000])
caxis([0 1e-2])

%% Pendekatan linear
tic
fprintf('Metode Pendekatan linear :\n')
Mall=M0;
yssdh=linear1d(M0,0,xobs);
Ermsaw=sqrt(1/numel(yobs)*sum((yssdh-yobs).^2));
for i=1:maxiter
    J(:,1)=(linear1d(M0-deltaM,0,xobs)-linear1d(M0,0,xobs))/deltaM;
    
    ysblm=linear1d(M0,0,xobs);
    M=M0+[J'*J+epsilon*eye(size(M0,1))]\J'*(ysblm-yobs);
    yssdh=linear1d(M,0,xobs);
    
    M0=M;
    Mall=cat(2,Mall,M0);
    Erms(i)=sqrt(1/numel(yobs)*sum((yssdh-yobs).^2));
    if Erms(i)<e
        fprintf('Number of iteration : %5i\n',i)
        fprintf('ERMS : %f\n',Erms(i))
        break
    end
end
waktu=toc;
if i==maxiter
    fprintf('Maximum iterration reached\n')
end
% Tampilan di console
fprintf('Hasil :\n')
fprintf('    a : %5.2f\n',M)
fprintf('Waktu running : %7.5fs\n\n',waktu)
% Plot ERMS per iterasi
figure
plot(Erms)
ylabel('ERMS')
xlabel('iterasi')
% Plot Spasial
ye=[Ermsaw Erms]+10;
figure
plot(A,t)
hold on
scatter(a,0,'*r')
ylabel('missfit'),xlabel('model parameter a')
axis([0 50 0 1000])
caxis([0 1e-2])
plot(Mall(:),ye,'--g')
scatter(Mall(1),ye(1),36,'k','filled','markeredgecolor','k')
for i=2:size(Mall,2)-1
    scatter(Mall(i),ye(i),36,'g','filled','markeredgecolor','k')
end
scatter(Mall(end),ye(end),36,'m','filled','markeredgecolor','k')
title('Pendekatan linear')

%% PSO
fprintf('Metode Particle Swarm Optimization :\n')
tic
% Initial Condition
Pbesta=ones(1,Nagent)*25;
Gbesta=ones(1);
va=zeros(1,Nagent);

% Initial fitness matrix
for i=1:Nagent
    ti(:,i)=linear1d(Apso(1,i),0,xobs);
    fitness(1,i)=sqrt(1/numel(yobs)*sum((ti(:,i)-yobs).^2));
end

% Initial global best
[rowg,colg,valg]=find(fitness==min(min(fitness)));
Gbesta=Apso(rowg(1),colg(1));

% Iterasi PSO
for k=1:15 
    for i=1:Nagent
        % Perhitungan posisi vektor A
        aa(k,i)=(Pbesta(i)-Apso(k,i))*c1*rand+(Gbesta-Apso(k,i))*c2*rand;
        va(k+1,i)=w*va(k,i)+aa(k,i);
        Apso(k+1,i)=Apso(k,i)+va(k+1,i);
        % Uji Error
        ti(:,i)=linear1d(Apso(k+1,i),0,xobs);
        fitness(k+1,i)=sqrt(1/numel(yobs)*sum((ti(:,i)-yobs).^2));
    end
    % Personal best
    [dum,Ib]=min(fitness,[],1);
    for i=1:Nagent
        Pbesta(i)=Apso(Ib(i),i);
    end
    % Global best
    [rowg,colg,valg]=find(fitness==min(min(fitness)));
    Gbesta=Apso(rowg(1),colg(1));
end
waktu=toc;
% Tampilan di console
fprintf('Global best :\n')
fprintf('    a : %5.2f\n',Gbesta)
fprintf('Waktu running : %7.5fs\n\n',waktu)
% Plot Kecepatan per iterasi
figure
plot(sqrt(va.^2))
xlabel('iterasi')
ylabel('kecepatan (km/iterasi)')
% Plot Spasial
for i=[1:16]
    figure
    plot(A,t)
    hold on
    scatter(a,0,'*r')
    ylabel('missfit'),xlabel('model parameter a')
    axis([0 50 0 1000])
    caxis([0 1e-2])
    scatter(Apso(i,:),fitness(i,:),200,'y','filled','markeredgecolor','k')
    text(Apso(i,:),fitness(i,:), num2str((1:Nagent)','%d'), 'horizontal','center')
    caxis([0 1e-2])
    saveas(gca,['iter' num2str(i) '.png'],'png')
end
%% Simulated Annealing
%model awal
tic
Msblm=rand*50;
MsaALLr=nan(1,2);
count=1;
for i=2:maxiter
    Mssdh=rand*50;
    
    ysblm=linear1d(Msblm,0,xobs);
    fnsblm=sqrt(1/numel(yobs)*sum((ysblm-yobs).^2));
    yssdh=linear1d(Mssdh,0,xobs);
    fnssdh=sqrt(1/numel(yobs)*sum((yssdh-yobs).^2));
    
    deltae=fnssdh-fnsblm;
    
    MsaALLb(i,:)=[Msblm fnsblm];
    
    %simple plot
    
%     plot(A,t)
%     hold on
%     scatter(a,0,'*r')
%     ylabel('missfit'),xlabel('model parameter a')
%     axis([0 50 0 1000])
%     caxis([0 1e-2])
%     scatter(Msblm,fnsblm+10,'ok')
    
    if deltae<=0
        Msblm=Mssdh;
        MsaALLr(count,:)=[Msblm,fnssdh];
        count=count+1;
    else
        P=exp(-deltae/T);
        R=rand;
        if R<P
            Msblm=Mssdh;
            MsaALLr(count,:)=[Msblm,fnssdh];
            count=count+1;
        end
    end
    
    T=T/log10(i);
end
waktu=toc;
% Tampilan di console
fprintf('Simulated Annealing :\n')
fprintf('    a : %5.2f\n',Msblm)
fprintf('Waktu running : %7.5fs\n\n',waktu)
% Plotting
figure
plot(A,t)
hold on
scatter(a,0,'*r')
ylabel('missfit'),xlabel('model parameter a')
axis([0 50 0 1000])
caxis([0 1e-2])
scatter(MsaALLb(:,1),MsaALLb(:,2),'ok','filled')
scatter(MsaALLr(:,1),MsaALLr(:,2),'or','filled')
    
%% Algoritma genetika
tic
fnga=nan(maxiterga,Npop);
AgaALL=nan(maxiterga,Npop);
fngas=nan(maxiterga,Npop);
AgaALLs=nan(maxiterga,Npop);
AgaALL(1,:)=round(Aga);

for i=2:maxiterga
    
    for j=1:Npop
        yga=linear1d(AgaALL(i-1,j),0,xobs);
        fnga(i-1,j)=sqrt(1/numel(yobs)*sum((yga-yobs).^2));
    end
    
    % Seleksi
    [Y,I]=sort(fnga(i-1,:));
    fngas(i-1,:)=fnga(i-1,I);
    AgaALLs(i-1,:)=AgaALL(i-1,I);
    
    
    % Reproduksi
    Agabin=dec2bin(AgaALLs(i-1,:),Nbin);
    Agabin(2*Popset+2:3*Popset+1,:)=Agabin(Popset+2:2*Popset+1,:);
    
    Agabinr(1,:)=Agabin(1,:);
    Agabinp=Agabin(randperm(size(Agabin,1)),:);
    for j=2:Npop
        Agabinr(j,:)=[Agabinp(j-1,1:co) Agabinp(j,co+1:end)];
    end
    
    % Mutasi
%     for j=randperm(Npop-1,2)+1
     for j=randi([1,Npop],1,2)
%         rm=randperm(Nbin,1);
        rm = randi([1,Nbin]);
        dump=str2double(Agabinr(j,rm));
        if dump==1
            dump=0;
        else
            dump=1;
        end
        Agabinr(j,rm)=num2str(dump);
    end
    AgaALL(i,:)=bin2dec(Agabinr)';
    
end
waktu=toc;
% Tampilan di console
fprintf('Algoritma Genetika :\n')
fprintf('    a : %5.2f\n',AgaALL(maxiterga,1))
fprintf('Waktu running : %7.5fs\n\n',waktu)
% Plot
figure
plot(fnga(:,1))
xlabel('iterasi')
ylabel('RMSE')
title('GA RMSE dari generasi terbaik')
figure
plot(AgaALL(:,1))
xlabel('iterasi')
ylabel('A')
title('GA generasi terbaik')
