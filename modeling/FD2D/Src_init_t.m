sigma = 0.5;  

gauss_pulse = exp(-t.^2/(2*sigma^2)) .* cos(2*pi*fc*t);

window = hann(length(t))';  

% 用窗函数做 taper
stf = gauss_pulse .* window;

if Src_disp==1
figure;
plot(t, stf);
xlabel('Time (s)');
ylabel('Amplitude');
title('Gaussian Pulse');

% 频谱
Fs = 1/dt;
N = length(stf);
f = Fs*(0:(N/2))/N;
Y = fft(stf);
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure;
plot(f, P1);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Gaussian Pulse');
xlim([0 60]);
grid on;
end

