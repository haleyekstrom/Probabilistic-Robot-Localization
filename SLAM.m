%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% SLAM %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Haley Ekstrom (81786834), June 2021 %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clc
clear all
delta_t=0.1;
t_f=307; %final time
k_f=t_f/0.1; %number of steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
study=1; %choose from {1:loop-closure study %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% or 2:monte carlo simulations} %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if study==1
    %loop-closure study process noise
    load w_x.mat;
    load w_y.mat;
    load w_theta.mat;
else
    %monte carlo simulations process noise
    w_x=normrnd(0,0.4*delta_t,[k_f,1]);
    w_y=normrnd(0,0.4*delta_t,[k_f,1]);
    w_theta=normrnd(0,0.05*delta_t,[k_f,1]);
end
%monte carlo parameters
M_run=15; %number of runs
total_number_of_noise_channels=5;
[indep_noise{1:total_number_of_noise_channels*M_run}] = RandStream.create('mrg32k3a','NumStreams',total_number_of_noise_channels*M_run);
for i_M=1:M_run
    noise_index_extract=(i_M-1)*total_number_of_noise_channels;
    for j=1:2
        z_noise{j,i_M} = randn(indep_noise{noise_index_extract+j},k_f,1); %measurement noise
    end
end
%landmark locations
nl=35;%number of landmarks
Xl{1}=[2;3];
Xl{2}=[-2;-21];
Xl{3}=[6;-3.5];
Xl{4}=[3;-8.5];
Xl{5}=[17;3];
Xl{6}=[18;-5];
Xl{7}=[52;-18];
Xl{8}=[25;-3.5];
Xl{9}=[32;14];
Xl{10}=[38;1];
Xl{11}=[36;9.5];
Xl{12}=[30;4];
Xl{13}=[52;24];
Xl{14}=[53;14];
Xl{15}=[63;21.5];
Xl{16}=[62;13];
Xl{17}=[75;22];
Xl{18}=[70;15];
Xl{19}=[82;21.5];
Xl{20}=[78;9];
Xl{21}=[42;35];
Xl{22}=[32;38];
Xl{23}=[47;27];
Xl{24}=[35;20];
Xl{25}=[45;40.5];
Xl{26}=[40;50];
Xl{27}=[60;2];
Xl{28}=[56;-6];
Xl{29}=[48;-7];
Xl{30}=[35;-12];
Xl{31}=[25;-15];
Xl{32}=[20;-19];
Xl{33}=[10;-14];
Xl{34}=[5;-16];
Xl{35}=[63;-11];
on_off=0; %controller on/off
for case_study=1:3 %measurement zone radii cases
    if case_study==1
        m_r=5;
        Col(case_study)='b';
    elseif case_study==2
        m_r=10;
        Col(case_study)='g';
    else
        m_r=15;
        Col(case_study)='r';
    end
    traj_name{case_study}=['estimated robot trajectory (r=',num2str(m_r),'m)'];
    landmark_name{case_study}=['r=',num2str(m_r),'m measurement zone'];
    for i_M = 1:M_run
        %initialization
        x_true{i_M}(1:3,1,case_study)= [0 0 0.01]'; %initial true robot states
        x_hat{i_M}(1:3,1,case_study)=[0,0,0.01]'; %initial predicted robot states
        P_hat{i_M}(:,:,1)=zeros(3+nl*2,3+nl*2); 
        P_hat{i_M}(1:3,1:3,1)=diag([1,1,(5*pi/180)^2]); %initial predicted robot error covariance
        Time(1)=0;
        j=1;
        for i = 1:35
            P_hat{i_M}(3+(2*i-1):3+2*i,3+(2*i-1):3+2*i,1)=100^2*eye(2,2);%initial predicted landmark error covariance
            x_hat{i_M}(3+(2*i-1):3+2*i,1,case_study)=[0;0]; %initial predicted landmark states
        end
        %noise statistics
        sigma_Q=diag([0.4*delta_t 0.4*delta_t 0.4*delta_t]);  %standard deviation of measurement error
        Q_robot=sigma_Q*sigma_Q; 
        Q=Q_robot; %process error covariance matrix
        for i = 1:35
            P_hat{i_M}(3+(2*i-1):3+2*i,3+(2*i-1):3+2*i,1)=100^2*eye(2,2);
            x_hat{i_M}(3+(2*i-1):3+2*i,1,case_study)=[0;0];
            Q=blkdiag(Q,zeros(2,2)); 
        end
        sigma_R=diag([0.1 (3*pi)/180]);%standard deviation of process error
        R=sigma_R*sigma_R; %measurement error covariance matrix
        for k=1:k_f
            %%%%%%%%%%%%%%%%%%%%%%%%%% simulator %%%%%%%%%%%%%%%%%%%%%%%%%%
            %control inputs
            if Time(k)<=20
                u(1,k)=1;
                u(2,k)=0;
            elseif (Time(k)>20 && Time(k)<=25)
                u(1,k)=1;
                u(2,k)=0.15;
            elseif (Time(k)>25 && Time(k)<=45)
                u(1,k)=1;
                u(2,k)=0;
            elseif (Time(k)>45 && Time(k)<=50)
                u(1,k)=1;
                u(2,k)=-0.15;
            elseif (Time(k)>50 && Time(k)<=75)
                u(1,k)=1;
                u(2,k)=0;
            elseif (Time(k)>75 && Time(k)<=90)
                u(1,k)=1;
                u(2,k)=-0.15;
            elseif (Time(k)>90 && Time(k)<=115)
                u(1,k)=1;
                u(2,k)=0;
            elseif (Time(k)>115 && Time(k)<=130)
                u(1,k)=1.2;
                u(2,k)=-0.15;
            elseif (Time(k)>130 && Time(k)<=160)
                u(1,k)=1.5;
                u(2,k)=0;
            elseif (Time(k)>160 && Time(k)<=190)
                u(1,k)=1.2;
                u(2,k)=-0.15;
            elseif (Time(k)>190 && Time(k)<=205)
                u(1,k)=1.2;
                u(2,k)=0.2;
            elseif (Time(k)>205 && Time(k)<=220)
                u(1,k)=1.5;
                u(2,k)=0;
            elseif (Time(k)>220 && Time(k)<=235)
                u(1,k)=1.5;
                u(2,k)=-0.2;
            elseif (Time(k)>235 && Time(k)<=285)
                u(1,k)=1.5;
                u(2,k)=0;
            elseif (Time(k)>285 && Time(k)<=290)
                u(1,k)=1.5;
                u(2,k)=-0.4;
            else
                u(1,k)=1;
                u(2,k)=0;
            end
            Time(k+1)=Time(k)+delta_t;
            for j=1:2
                v(j,k+1)=sigma_R(j,j)*z_noise{j,i_M}(k); %measurement noise
            end
            %%%%%%%%%%%%%%%%%%%%% extended kalman filter %%%%%%%%%%%%%%%%%%%%%%%
            %true robot states
            x_true{i_M}(1,k+1,case_study)=x_true{i_M}(1,k,case_study)+delta_t*u(1,k)*cos(x_true{i_M}(3,k,case_study))+w_x(k);
            x_true{i_M}(2,k+1,case_study)=x_true{i_M}(2,k,case_study)+delta_t*u(1,k)*sin(x_true{i_M}(3,k,case_study))+w_y(k);
            x_true{i_M}(3,k+1,case_study)=x_true{i_M}(3,k,case_study)+delta_t*u(2,k)+w_theta(k);
            for i_lm=1:nl
                x_true{i_M}(3+(2*i_lm-1):3+2*i_lm,k+1,case_study)=Xl{i_lm};%the landmarks are stationary 
            end
            %linearized system matrix 
            F=[];
            F=[1    0       -sin(x_hat{i_M}(3,k))*delta_t*u(1,k); %F_robot
               0    1        cos(x_hat{i_M}(3,k))*delta_t*u(1,k);
               0    0           1];
            %%%%%%%%%%%%%%%%%%%%%% prediction/propagation %%%%%%%%%%%%%%%%%%
            %robot state predictions
            x_bar{i_M}(1,k+1)=x_hat{i_M}(1,k,case_study)+delta_t*u(1,k)*cos(x_hat{i_M}(3,k,case_study));
            x_bar{i_M}(2,k+1)=x_hat{i_M}(2,k,case_study)+delta_t*u(1,k)*sin(x_hat{i_M}(3,k,case_study));
            x_bar{i_M}(3,k+1)=x_hat{i_M}(3,k,case_study)+delta_t*u(2,k);   
            for i_lm=1:nl
                x_bar{i_M}(3+(2*i_lm-1):3+2*i_lm,k+1,case_study)=x_hat{i_M}(3+(2*i_lm-1):3+2*i_lm,k,case_study); %landmark location predictions
                F=blkdiag(F,eye(2,2)); %F_landmark
            end
            P_bar{i_M}(:,:,k+1)=F*P_hat{i_M}(:,:,k)*F'+Q; %predicted error covariance
            %detecting landmarks in the measurment zone
            for i_l=1:nl
                range_to_landmark_i=sqrt((x_true{i_M}(1,k+1,case_study)-Xl{i_l}(1))^2+(x_true{i_M}(2,k+1,case_study)-Xl{i_l}(2))^2);
                if range_to_landmark_i<=m_r
                    landmark_detection_indicator{case_study,i_l}(k+1)=1;
                else
                    landmark_detection_indicator{case_study,i_l}(k+1)=0;
                end
            end
            z = [];
            z_hat = [];
            H = [];
            RR=[];
            detection_counter = 0;
            for i_l = 1:35
                range_to_landmark_i=sqrt((x_true{i_M}(1,k+1,case_study)-Xl{i_l}(1))^2+(x_true{i_M}(2,k+1,case_study)-Xl{i_l}(2))^2);
                if range_to_landmark_i<=m_r
                    detection_counter = detection_counter+1;
                    %true measurements
                    range_to_landmark=sqrt((x_true{i_M}(1,k+1)-x_bar{i_M}(3+(2*i_l-1),k+1)).^2+(x_true{i_M}(2,k+1)-x_bar{i_M}(3+(2*i_l),k+1)).^2); %true range measurement 
                    bearing_to_landmark=atan2(x_true{i_M}(2,k+1)-x_bar{i_M}(3+(2*i_l),k+1),x_true{i_M}(1,k+1)-x_bar{i_M}(3+(2*i_l-1),k+1))-x_true{i_M}(3,k+1); %true bearing measurement
                    if study==1
                        z_current_landmark=[range_to_landmark;bearing_to_landmark]+randn(2,1); %true measurement vector at current landmark
                    else
                        z_current_landmark=[range_to_landmark;bearing_to_landmark]+v(:,k+1); %true measurement vector at current landmark
                    end
                        %predicted measurements
                    range_pred=sqrt((x_bar{i_M}(1,k+1)-x_bar{i_M}(3+(2*i_l-1),k+1)).^2+(x_bar{i_M}(2,k+1)-x_bar{i_M}(3+2*i_l,k+1)).^2); %predicted range measurement
                    bearing_pred=atan2(x_bar{i_M}(2,k+1)-x_bar{i_M}(3+2*i_l,k+1),x_bar{i_M}(1,k+1)-x_bar{i_M}(3+(2*i_l-1),k+1))-x_bar{i_M}(3,k+1); %predicted bearing measurement    
                    z_estimated_to_current=[range_pred;bearing_pred]; %predicted measurement vector at current landmark
                    %linearized measurement matrices at current landmark
                    a = x_bar{i_M}(1,k+1)-x_bar{i_M}(3+(2*i_l-1),k+1);
                    b = x_bar{i_M}(2,k+1)-x_bar{i_M}(3+2*i_l,k+1);
                    H_current_robot=[((a)/(sqrt(a.^2+b.^2))),  ((b)/(sqrt(a.^2+b.^2))),  0; %linearized robot measurment matrix at current landmark
                        ((-b)/((1+(b/a).^2)*(a^2))),  (1/((1+(b/a).^2)*a)),  -1]; 
                    H_current_landmark=[((a)/(sqrt(a.^2+b.^2))), ((b)/(sqrt(a.^2+b.^2))),; %linearized landmark location measurement matrix at current landmark
                        ((b)/((1+(b/a).^2)*(a.^2))),  (-1/((1+(b/a).^2)*a))];
                    z=[z;z_current_landmark]; %true measurement vector
                    z_hat=[z_hat;z_estimated_to_current]; %predicted measurement vector
                    %linearized measurement matrix of current landmark
                    H_current=H_current_robot; 
                    for j=1:35
                        if j~=i_l
                            H_current=[H_current,zeros(2,2)];
                        else
                            H_current=[H_current,H_current_landmark];
                        end
                    end
                    H=[H;H_current]; %linearized measurement matrix
                    RR=blkdiag(RR,R); %process noise covariance matrix
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%% update %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if detection_counter>0
                S=H*P_bar{i_M}(:,:,k+1)*H'+RR; %innovation covariance
                K=P_bar{i_M}(:,:,k+1)*H'*inv(S); %kalman filter gain
                x_hat{i_M}(:,k+1,case_study)=x_bar{i_M}(:,k+1)+K*(z-z_hat); %corrected/updated state estimate
                P_hat{i_M}(:,:,k+1)=P_bar{i_M}(:,:,k+1)-K*S*K'; %corrected/updated error covariance
            else %no landmarks detected therefore no update/correction
                x_hat{i_M}(:,k+1,case_study)=x_bar{i_M}(:,k+1); %corrected/updated state estimate
                P_hat{i_M}(:,:,k+1)=P_bar{i_M}(:,:,k+1); %corrected/updated error covariance
            end
        end
    end
    loop_closure_error{case_study}=x_true{1}(1:2,1,case_study)-x_hat{1}(1:2,end,case_study); %loop closure error

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RMSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if study==2
    %computing RMSE of each state
        for i=1:3
            for k=1:k_f+1
                temp(i,k)=0;
                for i_M=1:M_run
                    error{i_M}(i,k)=(x_true{i_M}(i,k,case_study)-x_hat{i_M}(i,k,case_study))^2;
                    temp(i,k)=temp(i,k)+error{i_M}(i,k);
                end
                RMSE(i,k)=sqrt(temp(i,k)/(M_run));
            end
        end
        %plotting RMSE and absolute error in each individual run
        figure(3)
        tit=[num2str(M_run),' Monte Carlo runs'];
        set(gcf,'Position',[0 0 1000 1000]) %defining the size of the figure window
        for i_plot=1:3
            subplot(2,2,i_plot)
            hold on
            plot(Time(1:3071),RMSE(i_plot,:),Col(case_study),'LineWidth',2,'MarkerSize',10)
            xlabel('time')
            y_temp=['RMSE_',num2str(i_plot)];
            ylabel(y_temp);
            grid on
            set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2)
            title(tit)
            lgd=legend('r=5m','r=10m','r=15m');
            title(lgd,'Measurement Zones')
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%% loop-closure study plots %%%%%%%%%%%%%%%%%%%%%%%%%%%
if study==1
    for case_study=1:3
        hold on
        %plot 1: true and estimated robot trajectories over time
        figure(1) 
        set(gcf,'Position',[0 0 750 750]) %figure window size
        pbaspect([1 1 1]) %aspect ratio
        set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2)
        %plot true trajectory of robot
        true_mark=plot(x_true{1}(1,1,case_study),x_true{1}(2,1,case_study),'x','MarkerSize',12,'Color','k'); %mark the start point with x
        hold on
        true_trajectory=plot(x_true{1}(1,:,case_study),x_true{1}(2,:,case_study),'LineWidth',3,'Color','k'); %mark the start point with x
        hold on
        %plot estimated trajectory of robot
        est_mark=plot(x_hat{1}(1,1,case_study),x_hat{1}(2,1,case_study),'x','MarkerSize',12,'Color',Col(case_study));
        hold on
        est_trajectory{case_study}=plot(x_hat{1}(1,:,case_study),x_hat{1}(2,:,case_study),'LineWidth',3,'Color',Col(case_study));
        hold on
        %plot landmarks
        for j=1:nl
            xli=Xl{j}(1);
            yli=Xl{j}(2);
            landmark{j}=plot(xli,yli,'k*','MarkerSize',12);
        end
        axis([-20 100 -40 80])
        xlabel('x')
        ylabel('y')
        grid on
        hold on
        title('SLAM: Robot Trajectories')
        %plot 2: detected landmarks over time
        figure(2) 
        hold on
        for i_l=1:nl
            detected_landmark{case_study}(i_l)=plot(Time(2:k_f+1),landmark_detection_indicator{case_study,i_l}(2:k_f+1)*i_l,'Color',Col(case_study));
            hold on
        end
        set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2)
        xlabel('time')
        ylabel('landmark')
        yticks([1:nl])
        grid on
        title('SLAM: Detected Landmarks')
        axis([0 Time(k_f+1) 0 nl])
    end
    %plot legends
    legend([true_trajectory,est_trajectory{1},est_trajectory{2},est_trajectory{3},landmark{1}],'true robot trajectory',traj_name{1},traj_name{2},traj_name{3},'landmark')
    legend([detected_landmark{1}(1),detected_landmark{2}(1),detected_landmark{3}(1)],landmark_name{1},landmark_name{2},landmark_name{3})
end
   
  