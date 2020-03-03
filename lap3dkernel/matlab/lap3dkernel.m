% benchmark dense applying 1/r kernel between N and M pts in 3D.
% Barnett 3/3/20

clear
M = 1e4;  % targs
N = 1e4;  % sources

q = randn(N,1);  % source strengths
y = rand(3,N);    % source locs
x = rand(3,M);    % targ locs
u = zeros(M,1);  % the answer: pot at targs
itest = randi(M);  % which output to test
fprintf('test 1/r kernel in 3D. N=%d, M=%d...\n',N,M);

tic
for i=1:M   % outer loop over targs
  u(i) = sum(q'./sqrt((x(1,i)-y(1,:)).^2+(x(2,i)-y(2,:)).^2+(x(3,i)-y(3,:)).^2));
end
t = toc;
fprintf('targ-outer:\tu(test)=%.16g\t t=%.3g s\t %.3g Gpair/s\n',u(itest),t,N*M/t/1e9)

tic
u = zeros(M,1);
for j=1:N   % outer loop over src (adds to all targ pots each iteration)
  u = u + q(j)./sqrt((x(1,:)-y(1,j)).^2+(x(2,:)-y(2,j)).^2+(x(3,:)-y(3,j)).^2)';
end
t = toc;
fprintf('src-outer:\tu(test)=%.16g\t t=%.3g s\t %.3g Gpair/s\n',u(itest),t,N*M/t/1e9)

tic
A = 1./sqrt((x(1,:)'-y(1,:)).^2+(x(2,:)'-y(2,:)).^2+(x(3,:)'-y(3,:)).^2);
u = A*q;
t = toc;
fprintf('dense-matvec:\tu(test)=%.16g\t t=%.3g s\t %.3g Gpair/s\n',u(itest),t,N*M/t/1e9)

% maxes out at 0.1 Gpair/s on i7, even though uses all 8 threads.
