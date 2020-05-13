clear all
close all

filename=['../resolved/AdvectionOutC_1328_0_0.csv'];
R=csvread(filename,1);

filename=['AdvectionOut_320_07.csv'];
U=csvread(filename,1);

R=R(:,2:11);
U=U(:,2:11);

for f=1:10
    subplot(2,5,f)
    plot(R(:,f))
    hold on
    plot(U(:,f))
end
legend('Resolved','Unresolved')
