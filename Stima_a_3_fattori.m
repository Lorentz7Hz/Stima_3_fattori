%% Load and limit signal
load ("./5Vpp.mat");

clearvars -except signal_acq;
close all;
clc;

%% Dati
% frequenza campionamento
fs = 10^5;
% vettori tempo e frequenza
t_vector = 0:1/fs:(length(signal_acq)-1)/fs;
f_vector = (0:length(signal_acq)-1) * fs / length(signal_acq);
% range iniziale frequenze
start_freq = 999.000;
end_freq = 1001;
step = 0.001;
min_step = 1e-9; % Step minimo per evitare loop infiniti

%% MAIN stima a 3 fattori
while (end_freq - start_freq) > min_step
    
    freq_range = start_freq:step:end_freq;
    min_error = inf;
    best_freq = 0;

    for f = freq_range
        Di = [cos(2*pi*f*t_vector)', sin(2*pi*f*t_vector)', ones(length(signal_acq), 1)];
        stima = Di\signal_acq;
        sig = Di*stima;
        error = sqrt(sum((sig - signal_acq).^2));
        
        if error < min_error
            min_error = error;
            best_freq = f;
        end
    end

    % Aggiornamento della frequenza rispetto al miglior valore che approssima signal_acq
    start_freq = best_freq - step;
    end_freq = best_freq + step;
    step = step / 10; % Reduce step size
end

%display frequenza stimata
disp(['Frequenza stimata: ', num2str(best_freq), ' Hz']);

% calcolo parametri stimati
A0_est = stima(1);
B0_est = stima(2);
C0_est = stima(3);
A_est= sqrt(A0_est^2+B0_est^2);
phi_est =atan2(B0_est, A0_est);
signal_est_3params = A0_est*cos(2*pi*f*t_vector) + B0_est*sin(2*pi*f*t_vector) + C0_est;

%% Plots

figure(1);
plot(t_vector,signal_acq,t_vector, signal_est_3params);
xlim([0 0.1]);
grid on;
legend("acquired","estimated 3 params");

figure(2);
plot((sig)-signal_acq);
grid on;
legend("errore");

figure(3);
plot(db(fft(sig-signal_acq)));
grid on;
legend("segnale residuo");