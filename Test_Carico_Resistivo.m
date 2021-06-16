clear all
close all

% 
video = VideoWriter('C:\Users\mariu\Desktop\test resistenza', 'MPEG-4');
open(video); 

epsylon0 = 8.85e-12;
mu0 = 4 * pi * 1e-7;
C = (epsylon0 * 2.3 * 2 * pi)/log(10/2.8); 
L = mu0 / (2 * pi) * log(10/2.8); 
Lunghezza = 3;   
v = 1/sqrt(L*C);
fmax = 1.3e9; 
dz = v / (fmax * 50);
N = ceil(Lunghezza / dz);

        
dz = Lunghezza / (N-1);


T = 50e-9; 
dt = (0.5 * dz / v) / 500;

R_load = 5e1;
R_load = sqrt(L/C);

t0 = T/8; 
deltaT = 300e-12; 
Vin = @(t) normpdf(t, t0, deltaT);


figure;
fplot(@(x) Vin(x), [0, T])
title("legge V_{in} nel tempo");
xlabel("tempo [s]");
ylabel("tensione [V]")
grid()

Z = linspace(0, Lunghezza, N);
V = zeros(1, N);
I = zeros(1, N);

figure;
plotta_tensioni_correnti(Z, V, I);

pause()

tensioni = [];
correnti = [];

i = 0;
for t = 0:dt:T
    i = i + 1;
    
    V(1) = Vin(t);
    I(end) = V(end) / R_load;
    [V, I] = my_step(V, I, L, C, dt, dz);

    
    if (mod(i, 10000) == 0)
        disp(t/T);
        plotta_tensioni_correnti(Z, V, I);
        drawnow();
        frame = getframe(gcf);
        writeVideo (video, frame);
    end
    
    
    stride = 100;
    if (mod(i, stride) == 0)
        tensioni(1, i / stride) = V(1);
        correnti(1, i / stride) = I(1);
        tensioni(2, i / stride) = V(floor(end/2));
        correnti(2, i / stride) = I(floor(end/2));
        tensioni(3, i / stride) = V(end);
        correnti(3, i / stride) = I(end);
    end
    
end
close(video);

function [V, I] = my_step(V, I, L, C, dt, dz)
    N = numel(V);
   dv_dz = diff(V)/dz;
   di_dz = diff(I)/dz;
        
    I(1:end-1) = I(1:end-1) - (dv_dz / L) * dt;
    V(2:end) = V(2:end) - (di_dz / C) * dt;
end

function [] = plotta_tensioni_correnti(xx, V, I)
xlabel("spazio [m]")
yyaxis left
plot(xx, V);
ylabel("tensione [V]")
yyaxis right
plot(xx, I);
ylabel("corrente [A]")
yyaxis left
%grid();
end