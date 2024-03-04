%% Load and limit signal
load ("./5Vpp.mat");

clearvars -except signal_acq;
close all;
clc;

%% Dati
% frequenza di campionamento
fs = 10^5;
% vettori tempo e frequenza
t_vector = 0:1/fs:(length(signal_acq)-1)/fs;
f_vector = (0:length(signal_acq)-1) * fs / length(signal_acq);
% range iniziale frequenze
N = length(signal_acq);
frequencies = (0:N-1)*(fs/N);
fft_signal_acq = fft(signal_acq);
magnitude_spectrum = abs(fft_signal_acq(1:N/2+1));
[~, dominant_index] = max(magnitude_spectrum);
f0 = frequencies(dominant_index); % valore iniziale di f0 basato sulla frequenza dominante
start_freq = f0-1;
end_freq = f0+1;
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

csvwrite('signal_acq.csv', signal_acq);
csvwrite('sig.csv', sig);
% si consignlia di svolgere i seguenti linee di codice su altre piattaforme
% poichè il risultato vari ain base alla versione di matlab utilizzata
% rss = sum(signal_acq - sig).^2;
% tss = sum(signal_acq - mean(sig)).^2;
% r2 = 1 - rss/tss;
% disp(['Info convertita: ', num2str(r2*100), '%']);

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

% figure(2);
% plot((sig)-signal_acq);
% grid on;
% legend("errore");
% 
% figure(3);
% plot(db(fft(sig-signal_acq)));
% grid on;
% legend("segnale residuo");


%% Probabilità di errore

% Calculate the error between the acquired signal and the estimated signal
error_signal = sig - signal_acq;

% Create probability plot
figure;
probplot(error_signal);
title('Probability Plot of Errors in Comparison to Gaussian');

% Additional plot for visualization
figure;
subplot(2, 1, 1);
plot(error_signal);
title('Error Signal');
xlabel('Sample');
ylabel('Error');

subplot(2, 1, 2);
plot(db(fft(error_signal)));
title('FFT of Error Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');