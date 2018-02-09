clear
clc
tic;

% initialization 
n1=1;
n2=1.37;
R=((n1-n2)/(n1+n2))^2;
miu_a=0.1;
miu_s=100;
miu_t=miu_a+miu_s;
g=0.9;
weight0=1;
weightLoss_record=zeros(101,1); % depth resolved weight loss recording

% launching the photon packets
for n = 1:100000
    n
    direct_cos=[0,0,1]'; % initialize the direction cosine
    position=[0,0,0]';   % initialize the position
    weight=weight0-R;    % initialize the photo weight after specular reflection
    while weight>0.001
        weight1=weight;
        xi_1=rand(1);
        s=-log(xi_1)/miu_t;
        position=position+direct_cos*s; %update position

        xi_2=rand(1);
        cos_theta=1/2/g*(1+g^2-((1-g^2)/(1-g+2*g*xi_2))^2);
        sin_theta=sqrt(1-cos_theta^2);
        xi_3=rand(1);
        psi=2*pi*xi_3;
        if direct_cos(3)==1
            direct_cos=[sin_theta*cos(psi),sin_theta*sin(psi),cos_theta]';
        else
            direct_cos1=direct_cos;
            direct_cos(1:2,1)=[cos_theta,-sin(psi)*sin_theta/(sqrt(1-direct_cos1(3)^2)),sin_theta*direct_cos1(1)*cos(psi)/(sqrt(1-direct_cos1(3)^2));sin(psi)*sin_theta/(sqrt(1-direct_cos1(3)^2)),cos_theta,sin_theta*direct_cos1(2)*cos(psi)/(sqrt(1-direct_cos1(3)^2))]*direct_cos1;
            direct_cos(3)=-sin_theta*(sqrt(1-direct_cos1(3)^2))*cos(psi)+direct_cos1(3)*cos_theta;  %update direction cosine
        end
        
        % hitting the boundary
        if position(3)<0
            direct_cos=[1,0,0;0,1,0;0,0,-1]*direct_cos;  % re-update the direction cosine if hitting the boundary
            position=[1,0,0;0,1,0;0,0,-1]*position;      % re-update the position if hitting the boundary
            alph_i=acos(abs(position(3)));
            if alph_i>asin(n1/n2)                        % tell whether it is a total reflcetion
                R1=1;
            else
                alph_t=asin(n2*sin(alph_i)/n1);
                R1=1/2*((sin(alph_i-alph_t))^2/(sin(alph_i+alph_t))^2+(tan(alph_i-alph_t))^2/(tan(alph_i+alph_t))^2);
            end
            weight=weight*R1;                            % the reflected weight
        end

        delta_weight=miu_a/miu_t*weight;  % weight loss due to absorbation
        weight=weight-delta_weight; % update weight
        delta_w=weight1-weight;
        
        
        % the grid line separation and number of grid elements in the z
        % direction are 0.05 cm and 101, respectively.
        if position(3)<5
            N=fix(position(3)/0.05);
            ratio_N=position(3)/0.05-N;
            weightLoss_record(N+1:N+2,1)=weightLoss_record(N+1:N+2,1)+delta_w*[ratio_N;1-ratio_N];
        end
        
        % Russian roulette, tell whether a photo is dead or recoverd
        if weight<0.001              
            xi_4=rand(1);
            if xi_4<0.1
                weight=10*weight;
            end
        end
    end
    % the result of 1000 photo packets
    if n==1000
        weightLoss_record1=weightLoss_record/n;
        %weightLossSum_record1=zeros(101,1);
        %for i=1:101
            %weightLossSum_record1(i)=sum(weightLoss_record1(1:i));
        %end
        %weight_record1=1-weightLossSum_record1;
        figure;
        z=0:0.05:5;
        %plot(z,weight_record1);
        plot(z,weightLoss_record1/miu_a);
        title('Monte Carlo simulation with 1000 photon packets');
        xlabel('z [cm]');
        ylabel('Weight [AE]');
    end
    
    % the result of 10000 photo packets
    if n==10000
        weightLoss_record2=weightLoss_record/n;
        %weightLossSum_record1=zeros(101,1);
        %for i=1:101
            %weightLossSum_record1(i)=sum(weightLoss_record1(1:i));
        %end
        %weight_record1=1-weightLossSum_record1;
        figure;
        z=0:0.05:5;
        %plot(z,weight_record1);
        plot(z,weightLoss_record2/miu_a);
        title('Monte Carlo simulation with 10000 photon packets');
        xlabel('z [cm]');
        ylabel('Weight [AE]');
    end
    
    % the result of 100000 photo packets
    if n==100000
        weightLoss_record3=weightLoss_record/n;
        %weightLossSum_record1=zeros(101,1);
        %for i=1:101
            %weightLossSum_record1(i)=sum(weightLoss_record1(1:i));
        %end
        %weight_record1=1-weightLossSum_record1;
        figure;
        z=0:0.05:5;
        %plot(z,weight_record1);
        plot(z,weightLoss_record3/miu_a);
        title('Monte Carlo simulation with 100000 photon packets');
        xlabel('z [cm]');
        ylabel('Weight [AE]');
    end
end

%save('hk2.mat')
toc;