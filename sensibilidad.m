%% Inicializacion
clc; clear; close all; 
clc; clear; close all; randn('seed', 1e6); rand('twister',1e6);

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

%% Estimar los parametros con ruido en las mediciones
[minE_test, posminE_test] = min(E_test);
N = DN(posminE_test);
posN = find(DN==N);
Xn = X(1:DN(posN),:);
Yn = Y(1:DN(posN),:);
C = C_N(posN);
epsilon = epsilon_N(posN);
sigma = sigma_N(posN);
[Beta,b] = msvr(Xn,Yn,C,epsilon,sigma);

%% Parametros Sistema de potencia
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

%% Agregar ruido a las señales y calcular errores relativos
SNRd = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
SNRd = 50
f = 50;

for k = 1:length(SNRd)
    long = size(v_test,1);
    mini = 5;    
    for cycle = 1:round(max(time_test)*f)
        if long > round(Fm*(1/f)+mini)
            mfin = round(Fm*(1/f)+mini);
        else
            mfin = long;
        end
        timex = time_test(mini:mfin);
        a1 = 2*f*trapz(timex,v_test(mini:mfin)'.*cos(2*pi*f*timex));
        b1 = 2*f*trapz(timex,v_test(mini:mfin)'.*sin(2*pi*f*timex));
        Cv1(cycle,1) = sqrt(a1^2+b1^2);
        %
        a1 = 2*f*trapz(timex,i_test(mini:mfin)'.*cos(2*pi*f*timex));
        b1 = 2*f*trapz(timex,i_test(mini:mfin)'.*sin(2*pi*f*timex));
        Ci1(cycle,1) = sqrt(a1^2+b1^2); 
        %
        mini = mfin+1;
    end
    L = length(time_test);
    %
    powfund_v = mean(Cv1)^2/2;
    varnoise_v = powfund_v/10^(SNRd(k)/10);
    sv = sqrt(varnoise_v);
    ruido_v = sv*wgn(L,1,0);
    v_ruido = v_test + ruido_v;
    SNRv = snr(v_ruido,ruido_v);
    %
    powfund_i = mean(Ci1)^2/2;
    varnoise_i = powfund_i/10^(SNRd(k)/10);
    si = sqrt(varnoise_i);
    ruido_i = si*wgn(L,1,0);
    i_ruido = i_test + ruido_i;
    SNRi = snr(i_ruido,ruido_i);
    %
    %Calcular la componente fundamental de voltaje de Thevenin
    time_sim = max_time_test;
    datos_reales.time = time_test';
    datos.time = time_test';
    %
    datos_reales.signals.values = [time_test',i_ruido,v_ruido];
    sim estimar_vThevenin
    %
    mini = 5;
    long = size(v_the,1);
    vthe_1 = ones(long,1);
    %
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
    %
    fundamental.time = t;
    fundamental.signals.values = [t,vthe_1,v_pcc];
    %
    %Estimar los parametros y simular el sistema de potencia
    x_test = mean(abs(spectrogram(v_ruido/1000,fix(Fm*10/50),fix(Fm*9/50),F,Fm)),2);
    K = kernelmatrix(Xn,x_test',sigma); 
    ypred = (Beta'*K + b)';
    k1 = ypred(1); k2 = ypred(2); k3 = ypred(3); s = ypred(4); beta = ypred(5); w = ypred(6); Rc = ypred(7); Lc = ypred(8);
    %
    datos.signals.values = [time_test',i_ruido,v_ruido]; 
    sim validation
    %
    mini = 5;
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
        mini = mfin+1;
    end
    RErmsv = (abs(RMSv_real-RMSv_sim)./RMSv_real)*100;
    RErmsi = (abs(RMSi_real-RMSi_sim)./RMSi_real)*100;
    Results(k,:) = [round(SNRv,2),round(SNRi,2),round(mean(RErmsv),2),round(mean(RErmsi),2)];
end
Results
%% Figura para mostrar en el articulo
R = [6.5200    6.5100   19.1900   11.8800
     6.9100    6.9000   18.2100   10.5800
     7.3100    7.3000   17.0600   10.1700
     7.7300    7.7100    8.4100   12.2500
     8.1500    8.1300    8.4900   11.8900
     8.5700    8.5600    7.3000   11.2300
     9.0100    8.9900    7.0100   10.7900
     9.8900    9.8800    7.0200   9.7700
     10.8000   10.7900   6.8700   8.8100
     11.7300   11.7100   6.5500   8.0000
     12.6700   12.6500   6.4500   7.3500
     13.6200   13.6000   6.2800   6.8100
     14.5800   14.5700   6.0100   6.4200
     15.5500   15.5300   5.7800   6.1100
     16.5300   16.5100   5.7100   5.8500 
     17.5100   17.4900   5.7100   5.6300
     18.4900   18.4700   5.7300   5.4500
     19.4800   19.4600   5.7400   5.3100
     20.4700   20.4500   5.7400   5.2000
     25.4400   25.4200   5.8300   4.8200
     30.4400   30.4200   5.8300   4.6300
     35.4300   35.4100   5.7800   4.5500
     40.4300   40.4100   5.7300   4.5100
     50.4300   50.4100   5.6700   4.4800];
 
figure('Units','centimeters','Position',[9 6 8.5 5])
semilogx(R(:,1),R(:,3),'*-.','MarkerSize',6)
hold on
semilogx(R(:,2),R(:,4),'*-.','MarkerSize',6)
grid on
legend('Voltage','Current')
set(gca,'TickLabelInterpreter','latex','XLim',[6.5 50.6],'YLim',[4 20],'FontSize',8,'FontWeight','bold','FontName','Times');
ylabel({'Porcentual error [%]'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
xlabel({'SND (dB)'},'FontUnits','points','interpreter','latex','FontSize',8,'FontWeight','bold','FontName','Times')
matlab2tikz('figure_sensitive.tikz', 'height', '\figureheight', 'width', '\figurewidth');
%print('figure_sensitive','-depsc')