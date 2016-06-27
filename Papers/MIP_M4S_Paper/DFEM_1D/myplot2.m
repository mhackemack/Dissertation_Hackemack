function pout = myplot2(fid,phi,porder,dx,col)

% figure(fid);

ncells=length(dx);
f1=reshape(phi ,porder+1,ncells);f1=f1';
x1=0.;
hold on; grid on
for i=1:ncells,
    x2=sum(dx(1:i));
    xx=linspace(x1,x2,porder+1);
    if(porder>=1)
        xx=linspace(x1,x2,porder+1);
        y1=f1(i,:);
    else
        xx=[x1 x2];
        y1=f1(i)*ones(2,1);
    end
    pout = plot(xx,y1,col,'LineWidth',1);
%     hold on; grid on
    x1=x2;
end



y1=min(min(f1)); if(y1>0),y1=0;end
y2=max(max(f1)); y2=y2*1.05;
x1=0;
x2=sum(dx);
if(~isempty(y1>1)),y1=max(y1);end
if(~isempty(y2>1)),y2=max(y2);end
axis([x1 x2 y1 y2])
