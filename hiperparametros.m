%% Inicializacion
clc; clear; close all; 
[stat,struc] = fileattrib;
PathCurrent = struc.Name;
FolderName = 'hiperparametros';
PathFolder = [PathCurrent '\' FolderName];
mkdir([PathCurrent '/hiperparametros']);

%% Cargar datos
load ./datos_entrenamiento/N.mat;
load ./datos_entrenamiento/X.mat;
load ./datos_entrenamiento/Y.mat;
load ./datos_entrenamiento/F.mat

load ./datos_entrenamiento/time_train;
load ./datos_entrenamiento/max_time_train.mat;
load ./datos_entrenamiento/i_train.mat;
load ./datos_entrenamiento/v_train.mat;
load ./datos_entrenamiento/x_train.mat;

load ./datos_entrenamiento/time_test;
load ./datos_entrenamiento/max_time_test.mat;
load ./datos_entrenamiento/i_test.mat;
load ./datos_entrenamiento/v_test.mat;
load ./datos_entrenamiento/x_test.mat;

load ./datos_entrenamiento/Fm.mat
load ./datos_entrenamiento/Fs.mat;

%% Algoritmo ED
f = 50;
DN = [10,20,30,40,50,60,70,80,90,100,150,200,250,300]; save([PathFolder '/DN'],'DN');

gen = 30; save([PathFolder '/gen'],'gen');
NP = 20;
Fmut = 0.8;
CR = 0.5; 
rango_C = [1 100]';
rango_epsilon = [1e-6 1e-3]';
rango_sigma = [1 50]'; 

for n = 1:length(DN)    
    Xn = X(1:DN(n),:);
    Yn = Y(1:DN(n),:);   
   
    Em = zeros(NP,1);
    H = zeros(gen,3);
    
    Ho = compute_Ho(rango_C,rango_epsilon,rango_sigma,NP);
    
    time_sim = max_time_train;
    datos.time = time_train;
    datos.signals.values = [time_train,i_train,v_train];    
    for k = 1:gen
        for m = 1:NP 
            C = Ho(m,1);
            epsilon = Ho(m,2);
            sigma = Ho(m,3);
            [Beta,b] = msvr(Xn,Yn,C,epsilon,sigma);
            K = kernelmatrix(Xn,x_train',sigma);
            ypred = (Beta'*K + b)';
            k1 = ypred(1); k2 = ypred(2); k3 = ypred(3); s = ypred(4); beta = ypred(5); w = ypred(6); Rc = ypred(7); Lc = ypred(8);            
            sim modelo_EAF  
            E_Ho = computeOF(vreal_sim,vs_sim,t,f,Fs);
            %
            Ht = compute_Ht(NP,m,Ho,Fmut,CR);
            C = Ht(1);
            epsilon = Ht(2);
            sigma = Ht(3);
            [Beta,b] = msvr(Xn,Yn,C,epsilon,sigma);
            K = kernelmatrix(Xn,x_train',sigma);
            ypred = (Beta'*K + b)';
            k1 = ypred(1); k2 = ypred(2); k3 = ypred(3); s = ypred(4); beta = ypred(5); w = ypred(6); Rc = ypred(7); Lc = ypred(8); 
            sim modelo_EAF         
            E_Ht = computeOF(vreal_sim,vs_sim,t,f,Fs);
            
            if E_Ht <= E_Ho
                Ho(m,:) = Ht;
                Em(m,1) = E_Ht;
            else
                Em(m,1) = E_Ho;
            end
        end
        [E_train(n,k), posEm_min] = min(Em);
        H(k,:) = Ho(posEm_min,:);
    end    
    save([PathFolder '/E_train'],'E_train');    
    C = H(gen,1); 
    epsilon = H(gen,2); 
    sigma = H(gen,3);   
    [Beta,b] = msvr(Xn,Yn,C,epsilon,sigma);
    %
    K = kernelmatrix(Xn,x_test',sigma);
    ypred = (Beta'*K + b)';
    k1 = ypred(1); k2 = ypred(2); k3 = ypred(3); s = ypred(4); beta = ypred(5); w = ypred(6); Rc = ypred(7); Lc = ypred(8);     
    time_sim = max_time_test;
    datos.time = time_test';    
    datos.signals.values = [time_test',i_test,v_test];
    sim modelo_EAF    
    E_test(n) = computeOF(vreal_sim,vs_sim,t,f,Fs);
    %
    figure(1);
    plot(DN(n),E_train(n,gen),'r*');
    xlabel('Training Set Size');
    ylabel('Etrain');
    title(sprintf('Training Set Size = %d Fitness=%9.9f',DN(n),E_train(n,gen)));
    grid on;
    hold on;  
    C_N(n) = C;
    epsilon_N(n) = epsilon;
    sigma_N(n) = sigma;
end

save([PathFolder '/C_N'],'C_N');
save([PathFolder '/epsilon_N'],'epsilon_N');
save([PathFolder '/sigma_N'],'sigma_N');
save([PathFolder '/E_test'],'E_test');

%% Validacion de E_train y E_test
[minE_test, posminE_test] = min(E_test);
muestras = DN(posminE_test);
posN = find(DN==muestras);
Xn = X(1:DN(posN),:);
Yn = Y(1:DN(posN),:);
%
C = C_N(posN);
epsilon = epsilon_N(posN);
sigma = sigma_N(posN);
[Beta,b] = msvr(Xn,Yn,C,epsilon,sigma);
%
K = kernelmatrix(Xn,x_train',sigma);
ypred = (Beta'*K + b)';
k1 = ypred(1); k2 = ypred(2); k3 = ypred(3); s = ypred(4); beta = ypred(5); w = ypred(6); Rc = ypred(7); Lc = ypred(8);            
time_sim = max_time_train;
datos.time = time_train;
datos.signals.values = [time_train,i_train,v_train];
sim modelo_EAF
Comparar_Etrain = [computeOF(vreal_sim,vs_sim,t,f,Fs) E_train(posN,gen)];
%
K = kernelmatrix(Xn,x_test',sigma);
ypred = (Beta'*K + b)';
k1 = ypred(1); k2 = ypred(2); k3 = ypred(3); s = ypred(4); beta = ypred(5); w = ypred(6); Rc = ypred(7); Lc = ypred(8);            
datos.time = time_test';
time_sim = max_time_test;
datos.signals.values = [time_test',i_test,v_test];
sim modelo_EAF
Comparar_Etest = [computeOF(vreal_sim,vs_sim,t,f,Fs) E_test(posN)];

%% Figuras para reportar
figure('Units','centimeters','Position',[9 6 8.5 5.5])
plot(DN(1:12),E_train(1:12,gen),'*-.','MarkerSize',6,'MarkerEdgeColor',[0,0.447,0.741]);
grid on
set(gca,'TickLabelInterpreter','latex','XLim',[DN(1) DN(12)],'YLim',[4 6],'FontSize',8,'FontWeight','bold');
xlabel({'Training set size'},'FontUnits','points','Interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
ylabel({'Values of the OF'},'FontUnits','points','Interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
matlab2tikz('figure_Etrain.tikz', 'height', '\figureheight', 'width', '\figurewidth');

figure('Units','centimeters','Position',[9 6 8.5 5.5])
plot([1:gen],E_train(posN,:),'.k','MarkerSize',6,'MarkerEdgeColor',[0,0,0]);
grid on
set(gca,'TickLabelInterpreter','latex','XLim',[1 gen],'FontSize',8,'FontWeight','bold');
%title(sprintf('Iteration=%d, Fitness=%9.9f',gen,E_train(posN,gen)),'Interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times');
xlabel({'Iteration'},'FontUnits','points','Interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
ylabel({'Fitness'},'FontUnits','points','Interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
matlab2tikz('figure_OF.tikz', 'height', '\figureheight', 'width', '\figurewidth');
%print('figure_OF.tikz','-depsc')