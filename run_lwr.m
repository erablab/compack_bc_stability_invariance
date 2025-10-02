clc
clear
close all

data = Examples.lwr;

figure(1), hold on, grid on
plot(data.Vvec)
plot(data.Hvec)
plot(data.delta)
legend('V', 'h', 'delta')

figure(2), hold on, grid on, ylim([0,1])
plot(data.u_a)
plot(data.u_b)
legend('u_a', 'u_b')