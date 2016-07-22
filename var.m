
function []=jitt_demo

T = 1000;
m = 10;
ep = 0;

num_jitter=50;
num_runs=50;
jitter_width=1;   %this is delta
u=rand(num_runs,1);

%%%%%%%

clear pval pvalr pval_int pvalr_int
for ccc=1:num_runs

   orig_syn=0; 
    
    % neuron 1
    prob = zeros(T,1);
    odd = (1:2:T-1);
    even = (2:2:T);
    prob(odd) = 1-ep;
    prob(even) = ep;
    pos = rand(T,1) <= prob;
    n1 = find(pos');
    
    % neuron 2
    p = randperm(T); % generate a T array from 1 to T
    n2 = p(1:m);      % choose first m as m-length
    
    % compute synchrony
    
    orig_syn=synch_compute( n1,n2 );
    orig_synb=orig_syn+.5*rand(1);   % randomized synchrony

    % [basic] jitter, and tabulate synchrony counts

        syn_surr=[]; syn_surrb=[];
        for k=1:num_jitter
           
            % jitter spikes
            n1_jitt=n1;
            n2_jitt=n2+ jitter_width*(randi(3,1,length(n2))-2);
            % compute synchrony
            s=synch_compute(n1_jitt,n2_jitt);

            syn_surr(k)=s;
            syn_surrb(k)=s+.5*rand(1);   % store synchrony for surrogate j

        end
        
    % [interval] jitter, and tabulate synchrony counts

        syn_surr_int=[]; syn_surrb_int=[];
        for k=1:num_jitter

            % interval jitter (interval length jitter_width*2) spikes for n1
            n1_jitt_int=n1;
            %max( n1-n1_jitt )
            n2_jitt=(jitter_width)*floor( n2/(jitter_width) ) + (jitter_width)*(randi( 2,1,length(n2) )-1);
            
            % compute synchrony
            s=synch_compute( n1_jitt_int,n2_jitt);

            syn_surr_int(k)=s;
            syn_surrb_int(k)=s+.5*rand(1);   % store synchrony for surrogate j

        end
        
   % compute pvalues
        
   % p(X,R) for basic jitter
   pval(ccc)=(1+sum( syn_surr>=orig_syn))/(num_jitter+1);
   % p_C(X,R) for basic jitter
   pvalr(ccc)=(1+sum( syn_surrb>=orig_synb))/(num_jitter+1);
   % p(X,R) for interval jitter test
   pval_int(ccc)=(1+sum( syn_surr_int>=orig_syn))/(num_jitter+1);
   % p_C(X,R) for randomized interval jitter test
   pvalr_int(ccc)=(1+sum( syn_surrb_int>=orig_synb))/(num_jitter+1);
        


    if mod(ccc,200)==0
        orig_syn;,ccc
        binw=.01;
        subplot(3,2,1)
        hold off, histogram(pval,0:binw:1,'Normalization','probability'), title('p(X,R) basic jitter')
        hold on, plot(0:.005:1,binw*ones( size(0:.005:1)),'r-.')  % draw line
        subplot(3,2,3)
        hold off, histogram(pvalr,0:binw:1,'Normalization','probability'), title('p_C(X,R) basic jitter')
        hold on, plot(0:.005:1,binw*ones( size(0:.005:1)),'r-.')  % draw line
        subplot(3,2,2)
        hold off, histogram(pval_int,0:binw:1,'Normalization','probability'), title('p(X,R) interval jitter')
        hold on, plot(0:.005:1,binw*ones( size(0:.005:1)),'r-.')  % draw line
        subplot(3,2,4)
        hold off, histogram(pvalr_int,0:binw:1,'Normalization','probability'), title('p_C(X,R) interval jitter')
        hold on, plot(0:.005:1,binw*ones( size(0:.005:1)),'r-.')  % draw line
       
        subplot(3,2,5)
        hold off, histogram(u(1:ccc),0:binw:1,'Normalization','probability'), title('Uniform')
        hold on, plot(0:.005:1,binw*ones( size(0:.005:1)),'r-.')  % draw line
        pause(.001)
        
    end
                
end

function synch= synch_compute(n1,n2);
T =1000;
m = 10;
n_2=[];

N1 = zeros(T,1);
N1(n1) = 1;  % generate 1,0,1,0 alternate spike trail(for even T)

if n2(1)==0;
   n2 = n2(2:end);
end

for k=1:length(n2)
if n2(k) > T;
  n_2 = n2(1:k-1);
 break
end
end

N2 = zeros(T,1);
N2(n_2) = 1;
synch = dot(N1,N2);
