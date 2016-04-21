%Ujian Akhir Semester Inversi
%----------%Nama : Adi Wijaya Suchiana-----------------%
%----------%NIM  : 22313029----------------------------%

% harga maksimum dari f(x)=x+10 sin(5x)+7cos(4x)+25
% Genetic Algorithm
%%%----------------------------------------------------%%%


%input
clear all
Npop = 10; %jumlah populasi
maxIter = 3; % populasi maksimum
wF = zeros(1,Npop);
popGenBaru =zeros(1,Npop);
%%%----------------------------------------------------%%%


dum=0;
while dum==0
strt =1;
while strt==1
mulai = input('Start Genetic Algortihm (y/n)?','s'); % ketik y pada matlab
if strcmpi(mulai,'y')
   strt=0; 
elseif strcmpi(mulai,'n')
    dum=1;break;
else
    fprintf('Wrong Input'\n)
    strt=1;
end
end
if dum==1
    break
end

tic
%forward modeling untuk bantu plot
x=-5:0.001:5;
fxPlot = x + 10.*sin(5.*x)+7.*cos(4.*x)+25;
[ePlot,I]=max(fxPlot);
fprintf('Max pada x : %f\n',x(I));
fprintf('Nilai f Max : %f\n',ePlot)
plot (x,fxPlot);
hold on
plot (x(I),ePlot,'or')  
title ('plot fx')
xlabel('x')
ylabel('fx')
%%
%%%----------------------------------------------------%%%




% Langkah 1: populasi awal 10 bilangan acak
% pop = randi([-5,5],1,Npop);
pop = (5-(-5)).*rand(1,Npop)+ (-5);    %random real number dalam interval [-5,5]
iter =1;
l=1;
%Langkah 2: menentukan 5 pasangan terbaik berdasarkan fungsi f(x)
fx = pop + 10.*sin(5.*pop)+7.*cos(4.*pop)+25;
[fmaks,Im] = max(fx);   %harga maksimum
popM = pop(Im);         %x yang membuat f maksimum

%proses seleksi
while iter <maxIter
clear Ip fS popS
[furut,Iu] = sort(fx);
popU = pop(Iu);
totalF = sum(furut);       %jumlah total F untuk pembobotan
for i = 1:Npop
    wF(i) =  sum(furut(1:i))/totalF;   %normalisasi
end                     %
j=1;
for i=1:Npop
    if wF(i)>rand           %roulette wheel-->wF besar berarti fungsi f nya besar sehingga peluang terpilihnya besar
        Ip(j)=i;            %x yang lolos seleksi berdasarkan harga
        j=j+1;              %maksimum masing masing
    end
end
fS = furut(Ip);             
popS = popU(Ip);
Ipasangan = randi([1,size(Ip,2)],5,2); % indeks 5 pasangan terbaik
%Langkah 3 : reproduksi
R=rand;
popGenBaru(1,1:(size(Ipasangan)))= R.*popS(Ipasangan(:,1))+(1-R).*popS(Ipasangan(:,2));
popGenBaru(1,size(Ipasangan)+1:Npop)= (1-R).*popS(Ipasangan(:,1))+R.*popS(Ipasangan(:,2));

%Langkah 4: mutasi
mut=randi([1,10]);
popGenBaru(mut)=(5-(-5)).*rand+ (-5);


fxBaru = popGenBaru + 10.*sin(5.*popGenBaru)+7.*cos(4.*popGenBaru)+25;
[fmaksB,Im] = max(fxBaru);   %harga maksimum
popMB = popGenBaru(Im);      %x agar f maksimum
if fmaksB>fmaks
    fmaks=fmaksB;
    popM=popMB;
    fiter(l)=fmaks;
    popIter(l)=popM;
    l=l+1;
end

%Langkah 5 iterasi sebanyak 50 kali
iter=iter+1;
end

fprintf('Max pada x : %f\n',popM);
fprintf('Nilai f Max : %f\n',fmaks)
toc
dum=1;
%Langkah 6 plotting

plot (popM,fmaks,'*b')
legend('fungsi fx', 'nilai max','nilai max hasil GA')

figure
plot (x,fxPlot);
hold on
plot (x(I),ePlot,'or')  
title ('plot fx')
xlabel('x')
ylabel('fx')
title('perkembangan iterasi')
for i=1:length(l)
plot(popIter,fiter,'sk');
text(popIter,fiter, num2str((i)','%d'), 'horizontal','center')
end

end