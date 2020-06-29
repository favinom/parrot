function A=moveNodesPerfect(A,be,fe)

thickness=1.0e-4;
q=(0.5+0.75)/2;
centers=[0.5 q 0.75];

P=A(:,1);
points=unique(sort(P));

el2=be/2;
el4=be/4;
el8=be/8;

numEl=(length(points)-4)/2;

x1=linspace(0,centers(1)-thickness,el2+1);
f1=linspace(centers(1)-thickness,centers(1)+thickness,fe+1);

x2=linspace(centers(1)+thickness,centers(2)-thickness,el8+1);
f2=linspace(centers(2)-thickness,centers(2)+thickness,fe+1);

x3=linspace(centers(2)+thickness,centers(3)-thickness,el8+1);
f3=linspace(centers(3)-thickness,centers(3)+thickness,fe+1);

x4=linspace(centers(3)+thickness,1,el4+1);

x=[x1 f1(2:end) x2(2:end) f2(2:end) x3(2:end) f3(2:end) x4(2:end)];

y=abs(x-0.7);
[~,ii]=min(y);
x(ii)=0.7;

for coord=1:2
    for i=1:length(P)
        %if (mod(i,100))
        %    disp(i)
        %end
        A(i,coord)=x(points==A(i,coord));
    end
end
