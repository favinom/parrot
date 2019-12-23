clear all
close all

belist=[80,160,320,640];
felist=[1,2,4,8];

as=0;
us=0;

for bei=1:length(belist)
    for fei=1:length(felist)
        res=belist(bei)+3*felist(fei);
        filename=['AdvectionOut_',num2str(res),'_',num2str(us),'_',num2str(us),'.csv'];
        data{bei,fei}=csvread(filename,1);
        data{bei,fei}=data{bei,fei}(:,2:11);
    end
end


for r=1:size(data,1)
    figure
    clear mystring
    for fig=1:10
        subplot(2,5,fig)
        for c=1:size(data,2)
            plot(data{r,c}(:,fig))
            if (c==1)
                hold on
            end
            mystring{c}=['fe=',num2str(felist(c))];
        end
    end
    legend(mystring)
    sgtitle(['be=',num2str(belist(r))])
    filename=(['be',num2str(belist(r)),'.eps']);
    print(filename,'-depsc')
end



for c=1:size(data,2)
    figure
    clear mystring
    for fig=1:10
        subplot(2,5,fig)
        for r=1:size(data,1)
            plot(data{r,c}(:,fig))
            if (r==1)
                hold on
            end
            mystring{r}=['be=',num2str(belist(r))];
        end
    end
    legend(mystring)
    sgtitle(['fe=',num2str(felist(c))])
    filename=(['fe',num2str(felist(c)),'.eps']);
    print(filename,'-depsc')
end


figure
clear mystring
for fig=1:10
    subplot(2,5,fig)
    for d=1:min([size(data,1),size(data,2)])
        plot(data{d,d}(:,fig))
        if (d==1)
            hold on
        end
        mystring{d}=['be=',num2str(belist(d)),' fe=',num2str(felist(d))];
    end
end
legend(mystring)
filename=(['febe.eps']);
print(filename,'-depsc')
