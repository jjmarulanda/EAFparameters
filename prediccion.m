%% Inicializacion
clc; clear; close all; 

%% Cargar datos
load ./datos_entrenamiento/X.mat;
load ./datos_entrenamiento/Y.mat;
load ./datos_entrenamiento/F.mat

load ./datos_entrenamiento/time_test;
load ./datos_entrenamiento/max_time_test.mat;
load ./datos_entrenamiento/i_test.mat;
load ./datos_entrenamiento/v_test.mat;
load ./datos_entrenamiento/x_test.mat;

load ./datos_entrenamiento/Fm.mat
load ./datos_entrenamiento/Fs.mat;

load ./hiperparametros/C_N.mat
load ./hiperparametros/sigma_N.mat
load ./hiperparametros/epsilon_N.mat

load ./hiperparametros/gen.mat
load ./hiperparametros/DN.mat
load ./hiperparametros/E_train.mat
load ./hiperparametros/E_test.mat

%% Resultados con los datos de minimo error de test
[minE_test, posminE_test] = min(E_test);
N = DN(posminE_test);
posN = find(DN==N);
Xn = X(1:DN(posN),:);
Yn = Y(1:DN(posN),:);
C = C_N(posN);
epsilon = epsilon_N(posN);
sigma = sigma_N(posN);
[Beta,b] = msvr(Xn,Yn,C,epsilon,sigma);
K = kernelmatrix(Xn,x_test',sigma);
ypred = (Beta'*K + b)';
k1 = ypred(1); k2 = ypred(2); k3 = ypred(3); s = ypred(4); beta = ypred(5); w = ypred(6); Rc = ypred(7); Lc = ypred(8);

%% Parametros del sistema de potencia
Fs = 4*Fm;
f = 50;

%Thevenin
Rthe = 1.9;
Xthe = 5*37.7;

%T1 Transformer
St1 = 42e6;
V1t1 = 115000;
V2t1 = 11000;
aT1 = V1t1/V2t1;
Zt1_base_1 = V1t1^2/(St1/3); 
Lt1_pu = 0.12;
Rt1_pu = 0.012;

%Impedance
Zsec = 1e-6;

%EAF Transformer
Seaf = 30e6;
V1eaf = 11000;
V2eaf = 760;
aTeaf = V1eaf/V2eaf;
Zeaf_base_1 = V1eaf^2/(Seaf/3); 
Leaf_pu = 0.1;
Reaf_pu = 0.01;

%% Calcular el voltaje de Thevenin y simular el sistema de potencia
time_sim = max_time_test;
datos_reales.time = time_test';
datos_reales.signals.values = [time_test',i_test,v_test];
sim estimar_vThevenin

mini = 1;
long = size(v_the,1);
vthe_1 = ones(long,1);

for n = 1:round(max(t)*f)
    if long > round(Fs*(1/f)+mini)
        mfin = round(Fs*(1/f)+mini);
    else
        mfin = long;
    end
    timex = t(mini:mfin);
    a1 = (2*f)*trapz(timex,v_the(mini:mfin).*cos(2*pi*f*timex));
    b1 = (2*f)*trapz(timex,v_the(mini:mfin).*sin(2*pi*f*timex));
    C1 = sqrt(a1^2+b1^2);
    if a1>0 && b1>0
        theta_1 = atan(a1/b1);
    else
        if a1>0 && b1<0
            theta_1 = pi-atan(a1/abs(b1));
        else
            if a1<0 && b1<0
                theta_1 = pi+atan(a1/b1);
            else
                if a1<0 && b1>0
                    theta_1 = -atan(abs(a1)/b1);
                else
                    theta_1 = atan(a1/b1);
                end
            end
        end
    end
    vthe_1(mini:mfin) = C1*sin(2*pi*f*timex+theta_1);
    mini = mfin+1;
end

fundamental.time = t;
fundamental.signals.values = [t,vthe_1,v_pcc];

datos.time = time_test';
datos.signals.values = [time_test',i_test,v_test];  
sim validation

%% Espectro voltaje y corriente
x_pred = mean(abs(spectrogram(v_sim(1:4:end)/1000,fix(Fm*10/50),fix(Fm*9/50),F,Fm)),2);
RExv = (abs(x_test-x_pred)./x_test)*100; %error relativo
CRExv = [F(10)' RExv(10)]

xi_real = mean(abs(spectrogram(i_real(1:4:end)/100000,fix(Fm*10/50),fix(Fm*9/50),F,Fm)),2);
xi_pred = mean(abs(spectrogram(i_sim(1:4:end)/100000,fix(Fm*10/50),fix(Fm*9/50),F,Fm)),2);
RExi = (abs(xi_real-xi_pred)./xi_real)*100; %error relativo
CRExi = [F(10)' RExi(10)]

%% Metricas de funcionamiento y errores porcentuales
mini = 5;
long = size(v_real,1);

for n = 1:round(max(t)*f)
    if long > round(Fs*(1/f)+mini)
        mfin = round(Fs*(1/f)+mini);
    else
        mfin = long;
    end
    timex = t(mini:mfin);
    RMSv_real(n,1) = sqrt(f*trapz(timex,v_real(mini:mfin).*v_real(mini:mfin)));
    RMSi_real(n,1) = sqrt(f*trapz(timex,i_real(mini:mfin).*i_real(mini:mfin)));
    %
    RMSv_sim(n,1) = sqrt(f*trapz(timex,v_sim(mini:mfin).*v_sim(mini:mfin)));
    RMSi_sim(n,1) = sqrt(f*trapz(timex,i_sim(mini:mfin).*i_sim(mini:mfin)));
    %
    a1_real = 2*f*trapz(timex,v_real(mini:mfin).*cos(2*pi*f*timex));
    b1_real = 2*f*trapz(timex,v_real(mini:mfin).*sin(2*pi*f*timex));
    C1_real = sqrt(a1_real^2+b1_real^2)/sqrt(2);
    THDv_real(n,1) = 100*sqrt(RMSv_real(n)^2-C1_real^2)/C1_real;
    %
    a1_real = 2*f*trapz(timex,i_real(mini:mfin).*cos(2*pi*f*timex));
    b1_real = 2*f*trapz(timex,i_real(mini:mfin).*sin(2*pi*f*timex));
    C1_real = sqrt(a1_real^2+b1_real^2)/sqrt(2);
    THDi_real(n,1) = 100*sqrt(RMSi_real(n)^2-C1_real^2)/C1_real;
    %
    a1_sim = (2*f)*trapz(timex,v_sim(mini:mfin).*cos(2*pi*f*timex));
    b1_sim = (2*f)*trapz(timex,v_sim(mini:mfin).*sin(2*pi*f*timex));
    C1_sim = sqrt(a1_sim^2+b1_sim^2)/sqrt(2);
    THDv_sim(n,1) = 100*sqrt(RMSv_sim(n)^2-C1_sim^2)/C1_sim;    
    %
    a1_sim = (2*f)*trapz(timex,i_sim(mini:mfin).*cos(2*pi*f*timex));
    b1_sim = (2*f)*trapz(timex,i_sim(mini:mfin).*sin(2*pi*f*timex));
    C1_sim = sqrt(a1_sim^2+b1_sim^2)/sqrt(2);
    THDi_sim(n,1) = 100*sqrt(RMSi_sim(n)^2-C1_sim^2)/C1_sim;    
    %
    mini = mfin+1;
end

CompararRMSv = [RMSv_real RMSv_sim];
RErmsv = (abs(RMSv_real-RMSv_sim)./RMSv_real)*100; %error relativo valor RMS voltaje
Comparar_meanRMSv = [mean(RMSv_real) mean(RMSv_sim) (abs(mean(RMSv_real)-mean(RMSv_sim))./mean(RMSv_real))*100] 

CompararRMSi = [RMSi_real RMSi_sim];
RErmsi = (abs(RMSi_real-RMSi_sim)./RMSi_real)*100; %error relativo valor RMS corriente
Comparar_meanRMSi = [mean(RMSi_real)/1000 mean(RMSi_sim)/1000 (abs(mean(RMSi_real)-mean(RMSi_sim))./mean(RMSi_real))*100]

CompararTHDv = [THDv_real THDv_sim];
REthdv = (abs(THDv_real-THDv_sim)./THDv_real)*100; %error relativo valor THD voltaje
Comparar_meanTHDv = [mean(THDv_real) mean(THDv_sim)]

CompararTHDi = [THDi_real THDi_sim];
REthdi = (abs(THDi_real-THDi_sim)./THDi_real)*100; %error relativo valor THD corriente
Comparar_meanTHDi = [mean(THDi_real) mean(THDi_sim)]

[Pst_real, s_real] = flicker_sim(v_pcc/1000, 4*Fm, 50); %medicion de flicker
[Pst_sim, s_sim] = flicker_sim(v_pcc_sim/1000, 4*Fm, 50);
CPst = [Pst_real Pst_sim];
Er_pst = (abs(Pst_real-Pst_sim)./Pst_real)*100

%% Parametros del modelo para reportar
k1 = round(k1*10000,2); 
k2 = round(k2*10,2); 
k3 = round(k3*15,2); 
s = round(s*100,2); 
beta = round(beta*100,2); 
w = round(w*1,2); 
Rc = Rc*0.001; 
Lc = Lc*0.00001; 

%% Figuras para mostrar
ini = round(Fs*20/f);
fin = length(t);
ciclos_test = round(max(t)*50);

figure('Units','centimeters','Position',[9 6 8.5 5])
plot(t(ini:4:fin),v_real(ini:4:fin)/1000,'LineWidth',1)
hold on
plot(t(ini:4:fin),v_sim(ini:4:fin)/1000,'LineWidth',0.8)
grid on
legend('real','sim')
set(gca,'TickLabelInterpreter','latex','XLim',[t(ini) t(fin)],'FontSize',8,'FontWeight','bold','FontName','Times');
ylabel({'Voltage~[kV]'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
xlabel({'Time~[s]'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
matlab2tikz('figure_voltages.tikz', 'height', '\figureheight', 'width', '\figurewidth');
%print('figure_voltages','-depsc')

figure('Units','centimeters','Position',[9 6 8.5 5])
plot(t(ini:4:fin),i_real(ini:4:fin)/1000,'LineWidth',1)
hold on
plot(t(ini:4:fin),i_sim(ini:4:fin)/1000,'LineWidth',0.8)
grid on
legend('real','sim')
set(gca,'TickLabelInterpreter','latex','XLim',[t(ini) t(fin)],'YLim',[-1.2e2 1.2e2],'FontSize',8,'FontWeight','bold','FontName','Times');
ylabel({'Current~[kA]'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
xlabel({'Time~[s]'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
matlab2tikz('figure_currents.tikz', 'height', '\figureheight', 'width', '\figurewidth');
%print('figure_currents','-depsc')

figure('Units','centimeters','Position',[9 6 8.5 5])
plot([1:ciclos_test],RMSv_real,'LineWidth',1)
hold on
plot([1:ciclos_test],RMSv_sim,'LineWidth',0.8)
grid on
legend('real','sim')
set(gca,'TickLabelInterpreter','latex','XLim',[1 ciclos_test],'YLim',[0 500],'FontSize',8,'FontWeight','bold','FontName','Times');
ylabel({'Voltage~[kV]'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
xlabel({'Cicles'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
matlab2tikz('figure_RMSv.tikz', 'height', '\figureheight', 'width', '\figurewidth');
%print('figure_RMSv','-depsc')

figure('Units','centimeters','Position',[9 6 8.5 5])
plot([1:ciclos_test],RMSi_real/1000,'LineWidth',1)
hold on
plot([1:ciclos_test],RMSi_sim/1000,'LineWidth',0.8)
grid on
legend('real','sim')
set(gca,'TickLabelInterpreter','latex','XLim',[1 ciclos_test],'YLim',[0.2e2 0.7e2],'FontSize',8,'FontWeight','bold','FontName','Times');
ylabel({'Current~[kA]'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
xlabel({'Cicles'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
matlab2tikz('figure_RMSi.tikz', 'height', '\figureheight', 'width', '\figurewidth');
%print('figure_RMSi','-depsc')

figure('Units','centimeters','Position',[9 6 8.5 5])
plot(F,x_test,'LineWidth',1)
hold on
plot(F,x_pred,'--','LineWidth',0.8)
grid on
set(gca,'TickLabelInterpreter','latex','XLim',[0 650],'YLim',[0 65],'FontSize',8);
ylabel({'STFT'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
xlabel({'F~[Hz]'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
legend({'xreal','xest'},'interpreter','latex','FontSize',7,'FontWeight','bold','FontName','Times')
matlab2tikz('figure_spectrumV.tikz', 'height', '\figureheight', 'width', '\figurewidth');
%print('figure_spectrumV','-depsc')

figure('Units','centimeters','Position',[9 6 8.5 5])
plot(F,xi_real,'LineWidth',1)
hold on
plot(F,xi_pred,'--','LineWidth',0.8)
grid on
set(gca,'TickLabelInterpreter','latex','XLim',[0 650],'YLim',[0 90],'FontSize',8);
ylabel({'STFT'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
xlabel({'F~[Hz]'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
legend({'xreal','xest'},'interpreter','latex','FontSize',7,'FontWeight','bold','FontName','Times')
matlab2tikz('figure_spectrumI.tikz', 'height', '\figureheight', 'width', '\figurewidth');
%print('figure_spectrumI','-depsc')

figure('Units','centimeters','Position',[9 6 8.5 5])
plot(v_real(ini:4:fin)/1000,i_real(ini:4:fin)/100000,'LineWidth',1)
hold on
plot(v_sim(ini:4:fin)/1000,i_sim(ini:4:fin)/100000,'LineWidth',0.8)
grid on
legend('real','sim')
set(gca,'TickLabelInterpreter','latex','XLim',[-1 1],'YLim',[-1 1],'FontSize',8,'FontWeight','bold','FontName','Times');
ylabel({'Current~[kA]'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
xlabel({'Voltage~[kV]'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
matlab2tikz('figure_VI_characteristic.tikz', 'height', '\figureheight', 'width', '\figurewidth');
%print('figure_VI_characteristic','-depsc')