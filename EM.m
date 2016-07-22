clc;clear;

beta=input('Enter beta : ');
lam=input('Enter lambda vector : ');
m=input('Enter number of missing values in the x vector: ');
ep=input('Enter Max error in the parameter estomate : ');
MaxIter=input('Max number of allowable iteration: ');

[m1 n1]=size(beta); [m2 n2]=size(lam); 
if m1~=1 || n1~=1
    error('beta must be a scalar');
elseif beta<=0
    error('beta must be positive');
end
n=length(lam);

if m>=n
    error('Too many missing values');
end
if min(lam)<=0
    error('lambda must be positive');
end
x=poissrnd(lam((m+1):n));
y=poissrnd(beta*lam);
lam0=zeros(n,1);lam1=zeros(n,1);
lam0((m+1):n)=x;
beta0=mean((y((m+1):n))./x);
lam0(1:m)=y(1:m)/beta0;
iter=0;
while 1
    iter=iter+1;
    beta1=sum(y)/(sum(lam0(1:m))+sum(x));
    for i=1:m
        lam1(i)=(lam0(i)+y(i))/(beta1+1);
    end
    for i=(m+1):n
        lam1(i)=(x(i-m)+y(i))/(beta1+1);
    end
    if norm([beta0;lam0]-[beta1;lam1])<ep 
        beta_MLE=beta1;lam_MLE=lam1; break;
    elseif iter>MaxIter
        beta_MLE=beta1;lam_MLE=lam1; 
        warning('The alhorithm did not converge');break;
    else 
        lam0=lam1; beta0=beta1;
    end
end
disp('   ')
disp('**********')
disp('Number of iterations needed')
disp(iter)
disp('MLE for beta is')
disp(beta_MLE)
disp('MLE for lambda is')
disp(lam_MLE)
