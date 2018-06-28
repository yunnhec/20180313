clear;clc;
x=[6 4 0 5 0 0];
P=271;Wn=29;invWn=243;N=6;m=1;
m=m*3;
Wtmp=invWn;
wt=mod(Wtmp^(N/m),P);
k=0;
p=0;q=2;r=4;
tmp1=x(q+1);
tmp2=x(r+1);
x(r+1)=x(p+1)+wt^(k+m/3)*tmp1+wt^(k+2*m/3)*tmp2;
x(r+1)=mod(x(r+1),P);
x(q+1)=x(p+1)+wt^(k+2*m/3)*tmp1+wt^(k+4*m/3)*tmp2;
x(q+1)=mod(x(q+1),P);
x(p+1)=x(p+1)+wt^(k)*tmp1+wt^(k)*tmp2;
x(p+1)=mod(x(p+1),P);
m=m*2;
wt=mod(Wtmp^(N/m),P);
for k=0:m/2
    for p=k:m:N
        if q>=N-1 | p>=N-1
            break;
        else
            q=p+m/2;
            tmp1=x(q+1);
            x(q+1)=x(p+1)+wt^(k+m/2)*tmp1;
            x(q+1)=mod(x(q+1),P);
            x(p+1)=x(p+1)+wt^(k)*tmp1;
            x(p+1)=mod(x(p+1),P);
        end
    end
end
y=x;
x=[3 1 0 2 0 0];
P=271;Wn=29;invWn=243;N=6;m=1;
m=m*3;
Wtmp=invWn;
wt=mod(Wtmp^(N/m),P);
k=0;
p=0;q=2;r=4;
tmp1=x(q+1);
tmp2=x(r+1);
x(r+1)=x(p+1)+wt^(k+m/3)*tmp1+wt^(k+2*m/3)*tmp2;
x(r+1)=mod(x(r+1),P);
x(q+1)=x(p+1)+wt^(k+2*m/3)*tmp1+wt^(k+4*m/3)*tmp2;
x(q+1)=mod(x(q+1),P);
x(p+1)=x(p+1)+wt^(k)*tmp1+wt^(k)*tmp2;
x(p+1)=mod(x(p+1),P);
m=m*2;
wt=mod(Wtmp^(N/m),P);
for k=0:m/2
    for p=k:m:N
        if q>=N-1 | p>=N-1
            break;
        else
            q=p+m/2;
            tmp1=x(q+1);
            x(q+1)=x(p+1)+wt^(k+m/2)*tmp1;
            x(q+1)=mod(x(q+1),P);
            x(p+1)=x(p+1)+wt^(k)*tmp1;
            x(p+1)=mod(x(p+1),P);
        end
    end
end
x;
y;
x_y=mod(x.*y,P);
% ifft
m=1;
m=m*3;
Wtmp=invWn;
wt=mod(Wtmp^(N/m),P);
k=0;
p=0;q=2;r=4;
tmp1=x_y(q+1);
tmp2=x_y(r+1);
x_y(r+1)=x_y(p+1)+wt^(k+m/3)*tmp1+wt^(k+2*m/3)*tmp2;
x_y(r+1)=mod(x_y(r+1),P);
x_y(q+1)=x_y(p+1)+wt^(k+2*m/3)*tmp1+wt^(k+4*m/3)*tmp2;
x_y(q+1)=mod(x_y(q+1),P);
x_y(p+1)=x_y(p+1)+wt^(k)*tmp1+wt^(k)*tmp2;
x_y(p+1)=mod(x_y(p+1),P);
m=m*2;
wt=mod(Wtmp^(N/m),P);
for k=0:m/2
    for p=k:m:N
        if q>=N-1 | p>=N-1
            break;
        else
            q=p+m/2;
            tmp1=x_y(q+1);
            x_y(q+1)=x_y(p+1)+wt^(k+m/2)*tmp1;
            x_y(q+1)=mod(x_y(q+1),P);
            x_y(p+1)=x_y(p+1)+wt^(k)*tmp1;
            x_y(p+1)=mod(x_y(p+1),P);
        end
    end
end
% x_y=mod(x_y/N,P);