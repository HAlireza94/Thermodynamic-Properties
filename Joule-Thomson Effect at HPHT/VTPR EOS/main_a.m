clc;
n=500;
nsim=100;
alpha=-0.8;
phi=0.9;
sigmavi=0.363;sigmaeta=sigmavi;
thetatrue=[alpha phi sigmaeta];
randn('state',2^30-1);
epsilon=randn(n,1);
eta=randn(n,1);
x=zeros(n,1);
y=zeros(n,1);
x0=alpha/(1-phi^2);
x(1,1)=alpha+phi*x0+sigmaeta*eta(1,1);
y(1,1)=exp(x(1,1)*0.5)*epsilon(1,1);
for i=2:n
x(i,1)=alpha+phi*x(i-1,1)+sigmaeta*eta(i,1);
y(i,1)=exp(x(i,1)*0.5)*epsilon(i,1);
end
%%
theta=zeros(nsim,3);
thetam=zeros(nsim,3);
thetamse=zeros(nsim,3);
theta(1,1)=-0.5;%alpha;
theta(1,2)=0.5;%phi;
theta(1,3)=0.1;%sigmavi;
thetam(1,:)=theta(1,:);
x_ti=zeros(nsim,n);
x_ti(1,:)=x';
for k=2:nsim
theta(k,:)=theta(k-1,:);
x_ti(k,:)=x_ti(k-1,:);
theta(k,1)=alphasim(n,y,x_ti(k,:)',x0,theta(k,:));
theta(k,2)=phisim(n,y,x_ti(k,:)',x0,theta(k,:));
theta(k,3)=sigmaetasim(n,y,x_ti(k,:)',x0,theta(k,:));
thetam(k,:)=thetam(k-1,:)+theta(k,:);
x_ti(k,1)=x_1_sim(n,y,x_ti(k,:)',theta(k,:),x0);
for j=2:n-1
x_ti(k,j)=x_t_sim(n,y,x_ti(k,:)',theta(k,:),j);
end
x_ti(k,n)=x_n_sim(n,y,x_ti(k,:)',theta(k,:));
if (100*floor(k/100)==k)
disp(strvcat(['Sim=',num2str(k),', t=',num2str(j)]));
end
end
for i=1:1:3
thetam(:,i)=thetam(:,i)./(1:1:nsim)';
thetamse(:,i)=(thetam(:,i)-thetatrue(1,i)*ones(nsim,1)).^2;
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mi=min(x_ti,[],1)';
ma=max(x_ti,[],1)';
me=mean(x_ti,1)';
figure(1);
fill([(1:1:n)'; (n:-1:1)'],[mi((1:1:n)',1);ma((n:-1:1)',1)],[0.8 0.8 0.8]);
hold on;
plot(me,'-r');
plot(x,'-.b');
hold off;
legend('MCMC raw output',' MCMC Mean','True');
xlabel('Time');
set(gca,'Fontsize',14);
figure(2);
plot(me,'-r');
hold on;
plot(h,'-.k');
hold off;
legend(' MCMC Mean',' True');
xlabel('Time');
set(gca,'Fontsize',14);
figure(3);
plot(x_ti(:,[floor(n/2)-1]),'-k');
hold on;
plot(x([floor(n/2)-1],1)*ones(nsim,1),'-r');
hold off;
legend(strvcat(['x_{',num2str(floor(n/2)-1),'}']),'True');
title('MCMC rappresentative draws');
xlabel('MCMC Iterations');
set(gca,'Fontsize',14);
figure(4);
plot(x_ti(:,[floor(n/2)+1]),'-k');
hold on;
plot(x([floor(n/2)+1],1)*ones(nsim,1),'-r');
hold off;
legend(strvcat(['x_{',num2str(floor(n/2)+1),'}']),'True');
title('MCMC rappresentative draws');
xlabel('MCMC Iterations');
set(gca,'Fontsize',14);
figure(5)
hist(theta(:,1),100);
xlabel('\alpha');
%xlim([0 0.1]);
set(gca,'Fontsize',14);
figure(6);
hist(theta(:,2),100);
xlabel('\phi');
%xlim([1 1.2]);
set(gca,'Fontsize',14);
figure(7);
hist(theta(:,3),100);

xlabel('\sigma_{\eta}');
%xlim([0.1 0.2]);
set(gca,'Fontsize',14);
figure(8);
plot((1:1:size(theta,1))',theta(:,1));
hold on;
plot((1:1:size(theta,1))',ones(size(theta,1),1)*alpha,':r');
hold off;
ylabel('\alpha^{(t)}');
xlabel('MCMC Iterations (t)');
set(gca,'Fontsize',14);
figure(9);
plot((1:1:size(theta,1))',theta(:,2));
hold on;
plot((1:1:size(theta,1))',ones(size(theta,1),1)*phi,':r');
hold off;
ylabel('\phi^{(t)}');
xlabel('MCMC Iterations (t)');
set(gca,'Fontsize',14);
figure(10);
plot((1:1:size(theta,1))',theta(:,3));
hold on;
plot((1:1:size(theta,1))',ones(size(theta,1),1)*sigmavi,':r');
hold off;
ylabel('\sigma_{\eta}^{(t)}');
xlabel('MCMC Iterations (t)');
set(gca,'Fontsize',14);
figure(11);
plot((1:1:size(theta,1))',thetam(:,1));
hold on;
plot((1:1:size(theta,1))',ones(size(theta,1),1)*alpha,':r');
hold off;
ylabel('1/t \Sigma_{k=1}^{t} \alpha^{(k)}');
xlabel('MCMC Iterations (t)');
set(gca,'Fontsize',14);
figure(12);
plot((1:1:size(theta,1))',thetam(:,2));
hold on;
plot((1:1:size(theta,1))',ones(size(theta,1),1)*phi,':r');
hold off;
ylabel('1/t \Sigma_{k=1}^{t} \phi^{(k)}');
xlabel('MCMC Iterations (t)');
set(gca,'Fontsize',14);
figure(13);
plot((1:1:size(theta,1))',thetam(:,3));
hold on;
plot((1:1:size(theta,1))',ones(size(theta,1),1)*sigmavi,':r');
hold off;
ylabel('1/t \Sigma_{k=1}^{t} \sigma_{\eta}^{(k)}');
xlabel('MCMC Iterations (t)');
set(gca,'Fontsize',14);

