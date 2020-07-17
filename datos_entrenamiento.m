clc; clear; close all; randn('seed', 1e6); rand('twister',1e6);
[stat,struc] = fileattrib;
PathCurrent = struc.Name;
FolderName = 'datos_entrenamiento';
PathFolder = [PathCurrent '\' FolderName];
mkdir([PathCurrent '/datos_entrenamiento']);

%% Inicializacion
D = dlmread ('horno_Acasa.txt');

%% Vector de tiempo, voltajes y corrientes
t = D(:,1);
va = D(:,2);
vb = D(:,3);
vc = D(:,4);
ia = D(:,5);
ib = D(:,6);
ic = D(:,7);
Fm = round((t(end)-t(end-1))^(-1)); save([PathFolder '/Fm'],'Fm');

%% Seleccion de los datos en una de las fases
time = t;
i = ia;
v = va;

%% Separacion de los datos
ciclos = 100;
long_ciclos = round(Fm*(ciclos/50));      % numero de muestras en los ciclos considerados
long_time = length(t);                    % numero de muestras del vector de tiempo
div = round(long_time/long_ciclos);       % particiones posibles del vector de tiempo
new_time = time(1:long_ciclos*div);       % definicion del vector de tiempo considerando las particiones exactas
new_i = i(1:long_ciclos*div);             % definicion del vector de corriente considerando las particiones exactas
new_v = v(1:long_ciclos*div);             % definicion del vector de voltaje considerando las particiones exactas
T = reshape(new_time,[long_ciclos,div]);  % arreglo en cada columna de las particiones del vector de tiempo    
I = reshape(new_i,[long_ciclos,div]);     % arreglo en cada columna de las particiones del vector de corriente
V = reshape(new_v,[long_ciclos,div]);     % arreglo en cada columna de las particiones del vector de voltaje

%% Seleccion de los ciclos de corriente y tension
time = T(:,1);   
i_real = I(:,1); 
v_real = V(:,1); 

%% Frecuencias para caracterizar
F = [5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 150 200 250 300 350 400 450 500 550 600 650 700]; save([PathFolder '/F'],'F');

%% Datos de entrenamiento y test
porc_train = 70;

ciclos_train = round(porc_train*ciclos/100);    
time_train = time(1:round(Fm*ciclos_train/50)); save([PathFolder '/time_train'],'time_train');
max_time_train = 1.39990234375; save([PathFolder '/max_time_train'],'max_time_train');
i_train = i_real(1:round(Fm*ciclos_train/50));  save([PathFolder '/i_train'],'i_train'); 
v_train = v_real(1:round(Fm*ciclos_train/50));  save([PathFolder '/v_train'],'v_train');
x_train = mean(abs(spectrogram(v_train/1000,fix(Fm*10/50),fix(Fm*9/50),F,Fm)),2); save([PathFolder '/x_train'],'x_train'); 

ciclos_test = ciclos-ciclos_train; 
time_test = 0:1/Fm:ciclos_test/50; save([PathFolder '/time_test'],'time_test'); 
max_time_test = max(time_test); save([PathFolder '/max_time_test'],'max_time_test');
i_test = i_real(round(Fm*ciclos_train/50)+1:end); save([PathFolder '/i_test'],'i_test');
v_test = v_real(round(Fm*ciclos_train/50)+1:end); save([PathFolder '/v_test'],'v_test');
x_test = mean(abs(spectrogram(v_test/1000,fix(Fm*10/50),fix(Fm*9/50),F,Fm)),2); save([PathFolder '/x_test'],'x_test');

%% Datos de entrenamiento {X,Y}
k1min = 0.3;      k1max = 0.9;    %factor=10000
k2min = 0.1;      k2max = 0.5;    %factor=10
k3min = 0.3;      k3max = 0.8;    %factor=15
smin = 0.17;      smax = 0.39;    %factor=100
betamin = 0.14;   betamax = 0.20; %factor=100
wmin = 0.1;       wmax = 0.5;     %factor=1
Rcmin = 0.1;      Rcmax = 0.5;    %factor=0.001
Lcmin = 0.5;      Lcmax = 1;      %factor=0.00001

N = 1000; save([PathFolder '/N'],'N'); 
ymin = [k1min,k2min,k3min,smin,betamin,wmin,Rcmin,Lcmin];
ymax = [k1max,k2max,k3max,smax,betamax,wmax,Rcmax,Lcmax];
Y = computeY(N,ymin,ymax); save([PathFolder '/Y'],'Y');

Fs = 4*Fm; save([PathFolder '/Fs'],'Fs'); 

time_sim = max_time_train; 
datos.time = time_train;

for n = 1:N
    k1 = Y(n,1);   
    k2 = Y(n,2);
    k3 = Y(n,3); 
    s = Y(n,4);
    beta = Y(n,5);
    w = Y(n,6);
    Rc = Y(n,7);
    Lc = Y(n,8); 
    datos.signals.values = [time_train,i_train,v_train];
    sim modelo_EAF 
    X(n,:) = mean(abs(spectrogram(vs_sim(1:Fs/Fm:end)/1000,fix(Fm*10/50),fix(Fm*9/50),F,Fm)),2);
end
save([PathFolder '/X'],'X');