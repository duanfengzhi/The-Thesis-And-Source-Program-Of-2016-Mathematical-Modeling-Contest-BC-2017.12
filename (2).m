syms a b T0 p H C x y g x0 F0 T0 h0 v p2  rgg gg gr lgg lgt M mgg mgt rgt  th1 th2 th3  tha  T01 T12 T23 T34  f Fl y0 z(x) x1  f2
M=1000;g=9.8;mgg=10;mgt=100;p2=1025;rgg=0.025;lgg=1;rgt=0.15;lgt=1;p=28.12;
%mq=1200;
v=36
% % v=24,h0=0.7488;  v=12; h0= 0.7347;v=36,h0=0.7698,th6= 8.0541',b= 18.3167
%h0=0.7347;
%先让mq从初值增大，循环结束时的mq是是条件成立的最小值；(条件：两个角度要符合条件)
%再让mq从另一个初值减小，循环结束时的mq是最大值，即求得mq的范围。(条件：因为重物球质量越大，那两个角度越小，限制重物球上限的就是浮标吃水深度，设置h0<=1)
%最小值：mq=1780 h0=0.9440 h=18.0950 最大值：mq=1970 h0=1 h=17.8551
istrue=0;
%for mq=2000:-10:1800
    %for h0=0.8:0.01:1.0
    
   for h=16:0.1:20
 for mq=1600:20:4000
    for h0=0.7:0.01:1.2

     syms th4 th5 th6 T45 T56 x2
    
F0=1.25*v^2 *(2-h0); %风力

Ffl = p2*g*pi*h0; %3.1556*10^4*h   %浮力
th1 = atan((1.25*v^2 *(2-h0))/(3.1556*10^4*h0 - M*g));  %thi 第i节钢管的夹角
th2 = atan(F0/(Ffl-M*g + (p2*g*pi*rgg^2*lgg - mgg*g)));
th3 = atan(F0/(Ffl-M*g + 2*(p2*g*pi*rgg^2*lgg - mgg*g)));
th4 = atan(F0/(Ffl-M*g + 3*(p2*g*pi*rgg^2*lgg - mgg*g)));
T45 = F0 / sin(th4);    
th5 = atan(F0/(Ffl-M*g + 4*(p2*g*pi*rgg^2*lgg - mgg*g) + p2*g*pi*rgt^2*lgt - mgt*g - mq*g)); %钢桶倾角
th7= atan(F0/(Ffl-M*g + 4*(p2*g*pi*rgg^2*lgg - mgg*g) + p2*g*pi*rgt^2*lgt - mgt*g - mq*g - 22.05*7)); %无用
T56 = F0 / sin(th5);

ap = F0/(p*g);

th6 = solve('T56*sin(th5-th6) = T45*sin(th6-th4) + mq*g*sin(th6)',th6); %钢桶实际倾角 力矩
th6=th6(1);
%th6 =eval(th6(1));


f(x)=ap*(cosh(x/ap) - 1); %悬链线方程

x1=asinh(tan(pi/2 - th5))*ap
%z(x1)=vpa(subs(diff(f(x)),x,x1));
 %x1 = solve('atan(sinh(x1/ap))==pi/2 - th5',x1); %求锚链上端坐标x1
 %x1=eval(x1(1))
 
%x2=double(solve(ap*(tan(pi/2 - th5)-sinh(x2/ap))-22.05,x2));
x2=solve('ap*(sinh(x1/ap)-sinh(x2/ap))==22.05',x2);%根据锚链长度求锚链下端坐标x2
x2=eval(x2(1))


if x2<0 %锚链松弛
    l = eval(f(x1))
    b=0
    L=22.05-ap*(sinh(x1/ap)) %着地长度
else
    l = eval(f(x1)-f(x2))
    b=sinh(x2/ap) %锚点与海床的夹角
end


f2 =abs(eval(l + lgg*(cos(th1)+cos(th2)+cos(th3)+cos(th4))+lgt*cos(th6(1)) + h0 - h)/h) %计算深度误差

th6 =eval(th6(1))

if f2<0.05 &&th5>0  && rad2deg(b)<=16 && rad2deg(th6)<=5
    
   h= l + lgg*(cos(th1)+cos(th2)+cos(th3)+cos(th4))+lgt*cos(th6) + h0 ;        %预测海深
   r=abs(x1) + lgg*(sin(th1)+sin(th2)+sin(th3)+sin(th4))+lgt*sin(th6) %浮动半径
   istrue=1
   
if x2<0
    
 fplot(@(x)0,[-L 0],'b'); 
 hold on
 fplot(@(x)ap*(cosh(x/ap) - 1),[0 x1],'r'); 
 
else
  fplot(@(x)ap*(cosh(x/ap) - 1),[0 x2],'g'); 
  hold on
   fplot(@(x)ap*(cosh(x/ap) - 1),[x2 x1],'r'); 
    
end

 hold on
 fplot(@(x)(tan(pi/2-th6)*(x-x1)+f(x1)),[x1 x1+sin(th6)],'g'); 
 hold on
 fplot(@(x)tan(pi/2-th4)*(x-x1-sin(th6))+f(x1)+cos(th6),[x1+sin(th6) x1+sin(th6)+sin(th4)],'r'); 
 hold on
 fplot(@(x)tan(pi/2-th3)*(x-x1-sin(th6)-sin(th4))+f(x1)+cos(th6)+cos(th4),[x1+sin(th6)+sin(th4) x1+sin(th6)+sin(th4)+sin(th3)],'g'); 
 hold on
 fplot(@(x)tan(pi/2-th2)*(x-x1-sin(th6)-sin(th4)-sin(th3))+f(x1)+cos(th6)+cos(th4)+cos(th3),[x1+sin(th6)+sin(th4)+sin(th3) x1+sin(th6)+sin(th4)+sin(th3)+sin(th2)],'r'); 
 hold on
 fplot(@(x)tan(pi/2-th1)*(x-x1-sin(th6)-sin(th4)-sin(th3)-sin(th2))+f(x1)+cos(th6)+cos(th4)+cos(th3)+cos(th2),[x1+sin(th6)+sin(th4)+sin(th3)+sin(th2) x1+sin(th6)+sin(th4)+sin(th3)+sin(th2)+sin(th1)],'g'); 
 hold on

   
   %break;
end

   end
   
    
   %if istrue==1
    %break;
%end

 end
   end
 axis([-10 30 0 20]) 
     hold off


