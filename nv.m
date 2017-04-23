%~~ unconstraint problem for minimization whith shadowing path loss model~~~
clear,close all;
fd = fopen('resdvhop.txt','a+');
fprintf(fd,'-----%s-----\n',date);
fprintf(fd,'NodeAmount\tBeaconAmount\tR\tavg_err\tavg_time\tT\tcancltrl\n')
fd1 = fopen('tmpdvhop.txt','a+');
fprintf(fd1,'-----%s-----\n',date);
fprintf(fd1,'NodeAmount\tBeaconAmount\tR\terror\tEtime\ttrial\ttotal_count\tcancltrl\n')

nodeValues = 300;
beaconValues = (5:5:30)';
% rValues = 15;
maxTrial = 20;
total_err = 0;


BorderLength=100;       %Area of the network

for nodeValuesi = 1:size(nodeValues,1)
    for beaconValuesi = 1:size(beaconValues,1)
%          for rValuesi = 1:size(rValues,1)
             total_err = 0;
             total_time = 0;
             T=0;          
             cancltrl=0;
             
             for trial =    1:maxTrial
                total_count = 0;
                NodeAmount=nodeValues(nodeValuesi);               % Total Number of sensor node
                BeaconAmount=ceil(NodeAmount*beaconValues(beaconValuesi)/100);              % Number of Anchor node
               
                inputfilename = strcat('NodeDistr_',num2str(NodeAmount),'_',...
                    num2str(BeaconAmount),'_',num2str(trial),'.mat');%'..\data\..\TotalNode\',
                load(inputfilename);
%~~~~~~~~~~~~~~~~~~~~~~~~~~Selection of random number of anchor and sensor nodes~~~
%                 r1=.3; % communication ranging error
%                 r2=R*r1;
                clear r Sxy Beacon UN
%                 r=rand(NodeAmount,NodeAmount)*r2;
                Sxy=[[1:NodeAmount];C];
                Beacon=[Sxy(2,1:BeaconAmount);Sxy(3,1:BeaconAmount)];   %Anchor nodes coordinate
                UN=[Sxy(2,(BeaconAmount+1):NodeAmount);Sxy(3,(BeaconAmount+1):NodeAmount)];   %Unknown nodes coordinates
        
        f=5000000000;
        c=300000000;
        n2=4;   % path loss exponent
        l = 40;     %Communication radius of sensor node
        del= -78.5;   %threshold of received power
        sig=2;      % standard deviation
        Ps=20;      % sending power
        DOI=.4;
        PL_d0=20*log10((4*pi*f)/c); % free space path loss
        Proba = zeros(1,l);
        for r1=1:l
        clear Proba(r1)
            for i=1:100
                k=1+(rand-0.5)*DOI;  
                PL_d(r1,i)=k*(PL_d0+10*n2*log10(r1)-normrnd(0,sig));
                Pr(r1,i)=Ps-PL_d(r1,i);
            end
            a=mean(Pr(r1,:));

            Proba(r1)=0.5*(erfc((del-a)/((2^0.5)*sig)));
            % r(r1)=Proba(r1)*R;
        end
        clear R
        for r1 = l:-1:1
            if Proba(r1)>=0.90
                R=r1; 
                break;
            end
        end
                %~~~~~~~~~~~~~~~~~~~~~~~~~~Distance and hop count between all nodes to all nodes ~~~
clear h
                for i=1:NodeAmount
                    for j=1:NodeAmount
                        Dall(i,j)=((Sxy(2,i)-Sxy(2,j))^2+(Sxy(3,i)-Sxy(3,j))^2)^0.5;  %Distance between all nodes to all nodes
                        if (Dall(i,j)<=R)&&(Dall(i,j)>0)
                            h(i,j)=1;
                        elseif i==j
                            h(i,j)=0;
                        else h(i,j)=inf;
                        end
                    end
                end
                
%~~~~~~~~~~~~~~~~~~~~~~~~~Minimum hop count between all nodes to all nodes~~~
                for k=1:NodeAmount
                    for i=1:NodeAmount
                        for j=1:NodeAmount
                            if h(i,k)+h(k,j)<h(i,j)  %min(h(i,j),h(i,k)+h(k,j))
                                h(i,j)=h(i,k)+h(k,j);
                            end
                        end
                    end
                end
      %###############
%~~~~~~~~~~~~~~~~~~~~~~~hop count between anchor nodes by my method

Ef = R/BeaconAmount;
            
  clear h1 D1
                h1=h(1:BeaconAmount,1:BeaconAmount);        % Number of hops between Anchor node to Anchor node
                D1=Dall(1:BeaconAmount,1:BeaconAmount);     % Euclidean distance between Anchor node to Anchor
  
      for i=1:BeaconAmount
                      for j=1:BeaconAmount
                           if D1(i,j)==0;
                              Hup(i,j)=0;
                           else
                              Hup(i,j)=ceil(D1(i,j)/(R-Ef)); %Hop count between anchor nodes
                           end
                       end
                  end          
                
%~~~~~~~~~~~~~~~~~~~~~~~~~~~hopsize for every anchor pair

                for i=1:BeaconAmount
                    for j=1:BeaconAmount
                        if D1(i,j)==0
                            Hupsize(i,j)=0;
                        else
                            Hupsize(i,j)=D1(i,j)/Hup(i,j);
                        end

                    end
                end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~hopsize for every anchor node
                avghupsize=0;
                for i=1:BeaconAmount 
                    hupsize(i,1)=(sum(Hupsize(i,:))/(BeaconAmount-1)); % Hopsize of ith anchor node
                    hupsize1(i,1)=((R-hupsize(i,1))^2/hupsize(i,1))-(R/hupsize(i,1));% Improvement in hopsize
                    hupsize2(i,1)=hupsize(i,1)+hupsize1(i,1);
                    avghupsize=avghupsize+hupsize2(i,1);
                end
                avghup=avghupsize/BeaconAmount;%average hopsize
%           Eff=(R/5)-3;      
tic

clear hop1
hop1=h(1:BeaconAmount,(BeaconAmount+1):NodeAmount);  % Minimum number of hops between Anchor node to Unknown node

%~~~~~~~~~~~~~~~~~~~~~~~~~~Distance estimation between Unknown node to Anchor node by my method~~~
                
                clear Distance
                Distance= zeros(BeaconAmount,UNAmount);
                for i=1:UNAmount
                    for j=1:BeaconAmount
                        [minhop,mini] = min(hop1(:,i));    % Index of Nearest Anchor node to the unknown node
                        hopnum=hop1(j,i);                  % Minimum number of hops between Anchor node to Unknown node 
                        Distance(j,i)=avghup*hopnum;   % Distance between Anchor node and unknown node
                    end
                end
                d=Distance;
%~~~~~~~construction of matrix A, B and estimation of unknown nodes coordinate using quadratic programing 

                 clear  b e a1 a2 a3 A B hl x
                % x = zeros(2,UNAmount);
                b = zeros(BeaconAmount-1,1);
                t=10^6;
                
                for index = 1:UNAmount
                    e=(1./hop1(1:BeaconAmount,index));
                    e = diag([0;0;e]);
                    a1=zeros(2,BeaconAmount-1);
                    for i=1:2
                        for j=1:(BeaconAmount-1)
                            a1(i,j)=(Beacon(i,BeaconAmount)-Beacon(i,j))/(d(BeaconAmount,index)^2);
                        end
                    end
                    a1 = a1';
                    a2 = -diag(d(1:BeaconAmount-1,index))/(d(BeaconAmount,index)^2);
                    a3 = ones(BeaconAmount-1,1)/d(BeaconAmount,index);
                    A = 2*[a1 a2 a3];

                    for j = 1:BeaconAmount-1
                        b(j) = -Beacon(1,j)^2+Beacon(1,BeaconAmount)^2-Beacon(2,j)^2+...
                            Beacon(2,BeaconAmount)^2+d(j,index)^2-d(BeaconAmount,index)^2;
                    end
                    B = b/d(BeaconAmount,index)^2;
                    iter=0;
                    hl = e+t*A'*A;
%                     hl = e+t*(A'*A)+eye(size(A,2));
                    x(1:BeaconAmount+2,index) =t*(hl\A'*B);
                end
                clear X
                X = x(1:2,:);
                Etime = toc;
%~~~~~~~~~~~~~~~~~~ estimation of error and Accuracy~~~
                clear error
                error = zeros(1,UNAmount); 
                for i=1:UNAmount
                     if isnan(norm(X(:,i)))||norm(X(:,i))==Inf
                           error(1,i) = 0;   %total count contains NaN nodes in a trial
                           total_count = total_count + 1;
                     else
                         error(1,i)=(((X(1,i)-UN(1,i))^2+(X(2,i)-UN(2,i))^2)^0.5);
                     end
                end
                if UNAmount==total_count
                    cancltrl=cancltrl+1;
                    error = 0;
                else
                    T = T+total_count; % total no of node which do not localized in 20 trial
                    error=sum(error)/((UNAmount-total_count)*R);
                    total_err=total_err+error;
                    total_time = total_time + Etime;
                end 
             end
                avg_err = total_err/(maxTrial-cancltrl);
                avg_time = total_time/(maxTrial-cancltrl);

                fprintf(fd,'%g\t%g\t%g\t%g\t%g\t%g\t%g\n',NodeAmount,BeaconAmount,R,...
                   avg_err,avg_time,T,cancltrl); 
                fclose(fd);
                fd = fopen('resdvhop.txt','a+');
          end
%      end
end
 fclose(fd);       