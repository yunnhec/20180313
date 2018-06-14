function s = multi_str(s1,s2)

N = length(s1)+length(s2);
v = zeros(1,N);
v1 = zeros(1,N);
v2 = zeros(1,N);
v1(1:length(s1)) = fliplr(double(s1)-48); % 把輸入數字拆成多項式、再倒續
v2(1:length(s2)) = fliplr(double(s2)-48); % ex: s1 = 1234, double(s1)-48 = 1 2 3 4, fliplr(double(s1)-48) = 4 3 2 1
for i=1:N
    for j=1:i
        v(i) = v(i)+v1(j)*v2(i-j+1);
    end
end
for i=1:N-1
    v(i+1) = v(i+1)+floor(v(i)/10);
    v(i) = mod(v(i),10);
end
v = fliplr(v);
for i=1:N
    if v(i) ~= 0
        break;
    end
end
s = char(v(i:end)+48);

% an example of 1234*5678
% A=[4 3 2 1 0 0 0 0];
% B=[8 7 6 5 0 0 0 0];
% af=fft(A);
% bf=fft(B);
% cf=af.*bf;
% C=ifft(cf);
% cmf =0;
% for i=1:length(C)
%     cmf=cmf+C(i)*10^(i-1);
% end
% cmf
% cm=1234*5678;
% error = cmf-cm