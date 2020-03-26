A=csvread('../unresolved/AdvectionOut_40_08.csv',1);

t=A(:,1);
A=A(:,2:11);

for f=1:10
subplot(2,5,f)
plot(t,A(:,f))
end

