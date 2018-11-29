clc
clear

x = (-1.5:0.001:0.5)';
XX = zeros(length(x),1);
for i=1:length(x)
    XX(i,1)=HHSS(-2*x(i)-1);
end
syms s;
fplot(heaviside(-2*s-1),[-1.5,0.5]);
figure
plot(x,XX);