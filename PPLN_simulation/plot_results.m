% plot_results.m
% Simple plotting utilities to examine saved data (R_vs_dli.mat)
function plot_results()
    S = load(fullfile('outputs','R_vs_dli.mat'));
    dli = S.dli_vec; R = S.Rvec;
    figure;
    subplot(2,1,1);
    plot(dli*1e3, R,'-','LineWidth',1.6); xlabel('\Delta l_i (mm)'); ylabel('R (norm)'); grid on;
    title('Interferogram');
    subplot(2,1,2);
    % show FFT envelope to get bandwidth
    Y = abs(fftshift(fft(R)));
    f = linspace(-0.5,0.5,length(Y))*(1/(dli(2)-dli(1)));
    plot(f, Y); xlabel('spatial frequency (1/m)'); ylabel('|FFT|'); grid on;
    title('FFT of interferogram (qualitative)');
end
