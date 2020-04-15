close all
clear all

for i=1:8
    v{i}=zeros(4,3);
end

v{1}(1,:)=[0.05 0.25 0.5];
v{2}(1,:)=[0.5 0.05 0.05];
v{3}(1,:)=[0.05 1 0.5];
v{4}(1,:)=[0.05 1 0.48];
v{5}(1,:)=[0.17 1.9 0.7];
v{6}(1,:)=[0.23 1.9 0.7];
v{7}(1,:)=[0.77 1.9 0.7];
v{8}(1,:)=[0.83 1.9 0.7];

v{1}(2,:)=[0.95 0.25 0.5];
v{2}(2,:)=[0.5 0.05 0.95];
v{3}(2,:)=[0.95 1 0.5];
v{4}(2,:)=[0.95 1.0 0.48];
v{5}(2,:)=[0.17 1.9 0.3];
v{6}(2,:)=[0.23 1.9 0.3];
v{7}(2,:)=[0.77 1.9 0.3];
v{8}(2,:)=[0.83 1.9 0.3];

v{1}(3,:)=[0.95 2 0.5];
v{2}(3,:)=[0.5 0.3 0.95];
v{3}(3,:)=[0.95 2.2 0.85];
v{4}(3,:)=[0.95 2.2 0.14];
v{5}(3,:)=[0.23 2.2 0.3];
v{6}(3,:)=[0.17 2.2 0.3];
v{7}(3,:)=[0.77 2.2 0.3];
v{8}(3,:)=[0.83 2.2 0.3];

v{1}(4,:)=[0.05 2 0.5];
v{2}(4,:)=[0.5 0.3 0.05];
v{3}(4,:)=[0.05 2.2 0.85];
v{4}(4,:)=[0.05 2.2 0.14];
v{5}(4,:)=[0.23 2.2 0.7];
v{6}(4,:)=[0.17 2.2 0.7];
v{7}(4,:)=[0.77 2.2 0.7];
v{8}(4,:)=[0.83 2.2 0.7];

% for i=1:8
%     t=[v{i}; v{i}(1,:)];
%     plot3(t(:,1),t(:,2),t(:,3))
%     hold on
% end
% 
% legend('0','1','2','3','4','5','6','7')
format long

for i=1:8
    a{i}=norm(v{i}(2,:)-v{i}(1,:));
    b{i}=norm(v{i}(4,:)-v{i}(1,:));
    a{i}*b{i}*0.01
end

angle{1}=[0 0 0];
angle{2}=[0 90 0];
angle{3}=[0 0 16.2602];
angle{4}=[0 0 -15.8192];
angle{5}=[78.6901 -90 90 ]; % 
angle{6}=[-78.6901 -90 -90 ];
angle{7}=[0 -90 0];
angle{8}=[0 -90 0];

format short

for i=1:8
    c{i}=0.25*sum(v{i});
end

fx=['fx_string = ''',num2str(c{1}(1))];
for i=2:8
    fx=[fx,',',num2str(c{i}(1))];
end
fx=[fx,''''];

fy=['fy_string = ''',num2str(c{1}(2))];
for i=2:8
    fy=[fy,',',num2str(c{i}(2))];
end
fy=[fy,''''];

fz=['fz_string = ''',num2str(c{1}(3))];
for i=2:8
    fz=[fz,',',num2str(c{i}(3))];
end
fz=[fz,''''];

fd1=['fd1_string = ''',num2str(a{1})];
for i=2:8
    fd1=[fd1,',',num2str(a{i})];
end
fd1=[fd1,''''];

fd2=['fd2_string = ''',num2str(b{1})];
for i=2:8
    fd2=[fd2,',',num2str(b{i})];
end
fd2=[fd2,''''];

fd3=['fd3_string = ''','0.01'];
for i=2:8
    fd3=[fd3,',','0.01'];
end
fd3=[fd3,''''];

fa1=['fa1_string = ''',num2str(angle{1}(1))];
for i=2:8
    fa1=[fa1,',',num2str(angle{i}(1))];
end
fa1=[fa1,''''];

fa2=['fa2_string = ''',num2str(angle{1}(2))];
for i=2:8
    fa2=[fa2,',',num2str(angle{i}(2))];
end
fa2=[fa2,''''];

fa3=['fa3_string = ''',num2str(angle{1}(3))];
for i=2:8
    fa3=[fa3,',',num2str(angle{i}(3))];
end
fa3=[fa3,''''];

disp(fx)
disp(fy)
disp(fz)
disp(fd1)
disp(fd2)
disp(fd3)
disp(fa1)
disp(fa2)
disp(fa3)

return

for frac=1:8
p{frac}=zeros(4,3);
counter=0;
for j=-1:2:1
    for i=-1:2:1
        counter=counter+1;
        p{frac}(counter,:)=[i*a{frac}/2 j*b{frac}/2 0.0];
    end
end
p{frac}=rotateMy(p{frac},angle{frac});
p{frac}=p{frac}+repmat(c{frac},[4 1]);
p{frac}([3,4],:)=p{frac}([4,3],:);

t=[v{frac};v{frac}(1,:)];
q=[p{frac};p{frac}(1,:)];

plot3(t(:,1),t(:,2),t(:,3))
hold on
plot3(q(:,1),q(:,2),q(:,3))
xlabel('x')
ylabel('y')
zlabel('z')

max(max(abs(v{frac}-p{frac})))

end

function out=rotateMy(in,angles)

    angles=angles/180*pi;

    R1=zeros(3);
    R2=zeros(3);
    R3=zeros(3);
    
    R1(1,1)=cos(angles(1));
    R1(1,2)=-sin(angles(1));
    R1(2,1)=sin(angles(1));
    R1(2,2)=cos(angles(1));
    R1(3,3)=1.0;
    
    R2(1,1)=cos(angles(2));
    R2(1,3)=-sin(angles(2));
    R2(2,2)=1.0;
    R2(3,1)=sin(angles(2));
    R2(3,3)=cos(angles(2));
    
    R3(1,1)=1.0;
    R3(2,2)=cos(angles(3));
    R3(2,3)=-sin(angles(3));
    R3(3,2)=sin(angles(3));
    R3(3,3)=cos(angles(3));

    R=R1*R2*R3
    
    out=(R1*R2*R3*(in'))';

    
end