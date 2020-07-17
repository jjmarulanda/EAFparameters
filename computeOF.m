function OF = computeOF(vreal,vsim,t,f,Fs)
mini = 5;
long = size(vreal,1);
for cycle = 1:round(max(t)*f)
    if long > round(Fs*(1/f)+mini)
        mfin = round(Fs*(1/f)+mini);
    else
        mfin = long;
    end
    timex = t(mini:mfin);
    RMSv_real = sqrt(f*trapz(timex,vreal(mini:mfin).*vreal(mini:mfin)));
    RMSv_sim = sqrt(f*trapz(timex,vsim(mini:mfin).*vsim(mini:mfin)));
    Evrms(cycle,1) = RMSv_real-RMSv_sim; 
    %
    mini = mfin+1;
end
OF = sqrt(sum(Evrms.^2,1))/length(Evrms);
end


