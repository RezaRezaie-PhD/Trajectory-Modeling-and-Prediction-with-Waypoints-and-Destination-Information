

% Predicting trajectories with waypoint and destination information

clear all
% close all
clc

%--------------------

N1 = 50;
N2 = N1 + 60;
N3 = N2 + 40;

T = 15;  %7

F = [1 T 0 0;...
    0 1 0 0;...
    0 0 1 T;...
    0 0 0 1];

qx = 0.01;
qy = 0.01;

Q = [qx*T^3/3 qx*T^2/2 0 0;...
    qx*T^2/2 qx*T 0 0;...
    0 0 qy*T^3/3 qy*T^2/2;...
    0 0 qy*T^2/2 qy*T];

% q_0 = 0.005;
% Q = q_0*eye(4);


sigma = [10 10];
R = diag(sigma.^2);

H = [1 0 0 0;...
    0 0 1 0];

%================================= 

XN0_m_r = [10000;80;5000;30];
XN0_Cov_r = [10000 400 0 0;...
             400 100 0 0;...
             0 0 10000 400;...
             0 0 400 100];
         
XN1_m_r = [90000;70;30000;50];
XN1_Cov_r = [10000 400 0 0;...
             400 100 0 0;...
             0 0 10000 400;...
             0 0 400 100];

XN2_m_r = [170000;60;170000;60];
XN2_Cov_r = [10000 400 0 0;...
             400 100 0 0;...
             0 0 10000 400;...
             0 0 400 100];

         
XN3_m_r = [250000;90;200000;30];
XN3_Cov_r = [10000 400 0 0;...
             400 100 0 0;...
             0 0 10000 400;...
             0 0 400 100];

%........................

%--------------------

XN0_m_t = [10500;60;5500;50];

XN0_Cov_t = [100000 0 0 0;...
             0 10000 0 0;...
             0 0 100000 0;...
             0 0 0 10000];
                  
%-----------------------
         
%----------------------


XN1_m_t = [90500;50;30500;70];
         
XN1_Cov_t = [100000 0 0 0;...
             0 10000 0 0;...
             0 0 100000 0;...
             0 0 0 10000];
%----------------------

%---------------------------

XN2_m_t = [170500;40;170500;80];
          
XN2_Cov_t = [100000 0 0 0;...
             0 10000 0 0;...
             0 0 100000 0;...
             0 0 0 10000];
                   
%----------------------------------

%----------------------------------

XN3_m_t = [250500;70;200500;50];
         
XN3_Cov_t = [100000 0 0 0;...
             0 10000 0 0;...
             0 0 100000 0;...
             0 0 0 10000];
%----------------------------------

%------------------------------


kf = 4;

iteration = 1000;

%..................

% PN1_x = zeros(1,N1);
% PN1_vx = zeros(1,N1);
% 
% PN2_x = zeros(1,N2);
% PN2_vx = zeros(1,N2);
% 6
% PN3_x = zeros(1,N3);
% PN3_vx = zeros(1,N3);

%-----------------

EEp = zeros(iteration,N3+1);
EEv = zeros(iteration,N3+1);
EEz = zeros(iteration,N3+1);

tic

for iter=1:iteration

   %===================================Realization
%     iter
    %
    Xr = zeros(4,N3+1);
    
    %-------------------
    
    Dxn = chol(XN0_Cov_r,'lower');
    Xr(:,1) = XN0_m_r + Dxn*randn(4,1);
    
    %-------------------
    
    C_N1_N0 = [8000 200 0 0;...
               200 70 0 0;...
               0 0 8000 200;...
               0 0 200 70];
      
    
%     C_N1_N0 = [8000 0 0 0;...
%         0 0 0 0;...
%         0 0 8000 0;...
%         0 0 0 0];
       
    
%         C_N1_N0 = zeros(4,4);

    %++++++++++++++++++++++++++++++
    %................
    
    XN1p_Cov_r = XN1_Cov_r - C_N1_N0/XN0_Cov_r*C_N1_N0';
    
    Dxn = chol(XN1p_Cov_r,'lower');
    XN1p_m_r = XN1_m_r + C_N1_N0/XN0_Cov_r*(Xr(:,1) - XN0_m_r);
    Xr(:,N1+1) = XN1p_m_r + Dxn*randn(4,1);
    
    %--------------------
    
    C_N2_N1 = [8000 200 0 0;...
               200 70 0 0;...
               0 0 8000 200;...
               0 0 200 70];
      
      
    
%     C_N2_N1 = [8000 0 0 0;...
%         0 0 0 0;...
%         0 0 8000 0;...
%         0 0 0 0];
    
    
%         C_N2_N1 = zeros(4,4);
    
    %................
    
    XN2p_Cov_r = XN2_Cov_r - C_N2_N1/XN1_Cov_r*C_N2_N1';
    
    Dxn = chol(XN2p_Cov_r,'lower');
    Xr(:,N2+1) = XN2_m_r + C_N2_N1/XN1_Cov_r*(Xr(:,N1+1) - XN1_m_r) + Dxn*randn(4,1);
    
    %----------------------------------
    
    C_N3_N2 = [8000 200 0 0;...
               200 70 0 0;...
               0 0 8000 200;...
               0 0 200 70];
      
      
    
    
%     C_N3_N2 = [8000 30 0 0;...
%                30 60 0 0;...
%                0 0 8000 30;...
%                0 0 30 60];
   
    
%         C_N3_N2 = zeros(4,4);
    
    %............................
    
    XN3p_Cov_r = XN3_Cov_r - C_N3_N2/XN2_Cov_r*C_N3_N2';
    
    Dxn = chol(XN3p_Cov_r,'lower');
    Xr(:,N3+1) = XN3_m_r + C_N3_N2/XN2_Cov_r*(Xr(:,N2+1) - XN2_m_r) + Dxn*randn(4,1);
    
    %++++++++++++++++++++++++++++++++++++++++++
    
    for k = 1:N3-1
        
        
        if k<N1
            
            %           k
            
            %-----------CN|k
            
            CNk = zeros(4,4);
            
            for ii=0:N1-k-1
                
                CNk = CNk + F^ii*Q*(F^ii)';
                
            end
            
            %============
            %----------
            
            Gk = Q - Q*F^(N1-k)'/(CNk + F^(N1-k)*Q*F^(N1-k)')*F^(N1-k)*Q;
            DG = chol(Gk,'lower');
            
            Gk_N = Gk*F^(N1-k)'/(CNk);
            Gk_km1 = F - Gk_N*F^(N1-k+1);
            
            %..........
            
            Xr(:,k+1) = Gk_km1*Xr(:,k) + Gk_N*Xr(:,N1+1) + DG*randn(4,1);
            
            %----------------
            
        elseif k>N1 && k<N2
            
            %-----------CN|k
            
            CNk = zeros(4,4);
            
            for ii=0:N2-k-1
                
                CNk = CNk + F^ii*Q*(F^ii)';
                
            end
            
            %============
            %----------
            
            Gk = Q - Q*F^(N2-k)'/(CNk + F^(N2-k)*Q*F^(N2-k)')*F^(N2-k)*Q;
            DG = chol(Gk,'lower');
            
            Gk_N = Gk*F^(N2-k)'/(CNk);
            Gk_km1 = F - Gk_N*F^(N2-k+1);
            
            %..........
            
            Xr(:,k+1) = Gk_km1*Xr(:,k) + Gk_N*Xr(:,N2+1) + DG*randn(4,1);
            
            %----------------
            %============
            
            
        elseif k>N2 && k<N3
            
            %-----------CN|k
            
            CNk = zeros(4,4);
            
            for ii=0:N3-k-1
                
                CNk = CNk + F^ii*Q*(F^ii)';
                
            end
            
            %============
            %----------
            
            Gk = Q - Q*F^(N3-k)'/(CNk + F^(N3-k)*Q*F^(N3-k)')*F^(N3-k)*Q;
            DG = chol(Gk,'lower');
            
            Gk_N = Gk*F^(N3-k)'/(CNk);
            Gk_km1 = F - Gk_N*F^(N3-k+1);
            
            %..........
            
            Xr(:,k+1) = Gk_km1*Xr(:,k) + Gk_N*Xr(:,N3+1) + DG*randn(4,1);
            
            %----------------
            %============
            
            
        end 
        
    end
    
    %===============================
%
%     %-------------------------------------------- Measurement
%
    Dxn = chol(R,'lower');
    Z = [Xr(1,2:end);Xr(3,2:end)] + Dxn*randn(2,N3);
%
%     %--------------------------------------------
%
%     %==========================================================TRACKING
%
%     %--------------------------------Initializaation

% %+++++++++++++++++++Mismatched Cross Covariance of Waypoints

% C_N1_N0 = [6000 100 0 0;...
%            100 50 0 0;...
%            0 0 6000 100;...
%            0 0 100 50];
%     
%        C_N2_N1 = [6000 100 0 0;...
%            100 50 0 0;...
%            0 0 6000 100;...
%            0 0 100 50];
%        
%        C_N3_N2 = [6000 100 0 0;...
%            100 50 0 0;...
%            0 0 6000 100;...
%            0 0 100 50];
% %        
%---------------------------------
%      C_N1_N0 = zeros(4,4);
%     C_N2_N1 = zeros(4,4);
%     C_N3_N2 = zeros(4,4);
    
    %------------------------------
    
%     C_N1_N0 = C_N1_N0/XN0_Cov_r*XN0_Cov_t;
%     C_N2_N1 = C_N2_N1/XN1_Cov_r*XN1_Cov_t;
%     C_N3_N2 = C_N3_N2/XN2_Cov_r*XN2_Cov_t;
% %     
%     C_N1_N0 = [70000 0 0 0;...
%                0 6000 0 0;...
%                0 0 70000 0;...
%                0 0 0 6000];
% %            
% 
%      C_N2_N1 = C_N1_N0;
%      C_N3_N2 = C_N1_N0;
    %------------------------------

    
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
    Ys = [XN0_m_t;XN1_m_t];
    PpN0 = XN0_Cov_t;
    PpN1 = XN1_Cov_t;
%     Pyp = blkdiag(Pp0,PpN);
    Pys = [PpN0 C_N1_N0';
           C_N1_N0 PpN1];

       
%     %--------------------
%
    Yh = zeros(8,N3);
    Xh = zeros(4,N3+1);

    Hy = H*[eye(4) zeros(4,4)];
%     %======================================= Main Loop

    Yh(:,1) = Ys;
    Xh(:,1) = XN0_m_t;
    
    for k = 1:N3


        if k<kf+1
            
            
            
            %---------------------------Prediction from k-1 to k
            
            %............................Dynamic model coefficients
            %............. CN|k
            
            CNk = zeros(4,4);
            
            for ii=0:N1-k-1
                
                CNk = CNk + F^ii*Q*(F^ii)';
                
            end
            %.............
            
            Gk = Q - Q*F^(N1-k)'/(CNk + F^(N1-k)*Q*F^(N1-k)')*F^(N1-k)*Q;
            
            Gk_N = Gk*F^(N1-k)'/(CNk);
            Gk_km1 = F - Gk_N*F^(N1-k+1);
            
            
            
            Gy1 = horzcat(Gk_km1,Gk_N);
            Gy2 = horzcat(zeros(4,4),eye(4));
            Gy = vertcat(Gy1,Gy2);
            
            GCy = blkdiag(Gk,zeros(4,4));
            
            %.............
            
            Yp = Gy*Ys;
            Pyp = Gy*Pys*Gy' + GCy;
            
            %--------------------------------------------------
            
            Pyh = Pyp - (Pyp*Hy')/(Hy*Pyp*Hy'+R)*(Pyp*Hy')';
            
            Yh(:,k+1) = Yp + (Pyp*Hy')/(Hy*Pyp*Hy'+R)*(Z(:,k) - Hy*Yp);
            
            Xh(:,k+1) = [eye(4) zeros(4,4)]*Yh(:,k+1);
            Ph = [eye(4) zeros(4,4)]*Pyh*[eye(4) zeros(4,4)]';
            
            %..................
            %.................. For the next time
            
            Ys = Yh(:,k+1);
            Pys = Pyh;
            
            Ykf = Yh(:,k+1);
            Pykf = Pyh;
            
            %+++++++++++++++++++++
            
            
            
        elseif k>kf && k<N1
           

           %---------------------------Prediction from k-1 to k
            
            %............................Dynamic model coefficients
            %............. CN|k
            
            CNk = zeros(4,4);
            
            for ii=0:N1-k-1
                
                CNk = CNk + F^ii*Q*(F^ii)';
                
            end
            %.............
            
            Gk = Q - Q*F^(N1-k)'/(CNk + F^(N1-k)*Q*F^(N1-k)')*F^(N1-k)*Q;
            
            Gk_N = Gk*F^(N1-k)'/(CNk);
            Gk_km1 = F - Gk_N*F^(N1-k+1);
            
            
            
            Gy1 = horzcat(Gk_km1,Gk_N);
            Gy2 = horzcat(zeros(4,4),eye(4));
            Gy = vertcat(Gy1,Gy2);
            
            GCy = blkdiag(Gk,zeros(4,4));
            
            %.............
            
            Yp = Gy*Ys;
            Pyp = Gy*Pys*Gy' + GCy;
            
            %--------------------------------------------------
            
            Pyh = Pyp;
            
            Yh(:,k+1) = Yp ;
            
            Xh(:,k+1) = [eye(4) zeros(4,4)]*Yh(:,k+1);
            Ph = [eye(4) zeros(4,4)]*Pyh*[eye(4) zeros(4,4)]';
            
            %..................
            %.................. For the next time
            
            Ys = Yh(:,k+1);
            Pys = Pyh;
            
            
            
            %+++++++++++++++++++++
            

            %--------------------------------------------------------

        elseif k == N1


            Xh(:,N1+1) = [zeros(4,4) eye(4)]*Ykf;
            Ph_N1 = [zeros(4,4) eye(4)]*Pykf*[zeros(4,4) eye(4)]';

        %------------------------- N2
        
            G_N2_N1 = C_N2_N1/XN1_Cov_t;
            Xh(:,N2+1) = XN2_m_t + G_N2_N1*(Xh(:,N1+1) - XN1_m_t);
            GN2 = XN2_Cov_t - C_N2_N1/XN1_Cov_t*C_N2_N1';
            
            Ph_N2 = GN2 + G_N2_N1*Ph_N1*G_N2_N1';
            
            %-------------------N1 N2
            
            P_N2_N1 = G_N2_N1*Ph_N1;
            %-------------------------- y 
            
            Ys = [Xh(:,N1+1);Xh(:,N2+1)];
            Pys = [Ph_N1 P_N2_N1';
                   P_N2_N1 Ph_N2];
            
               %-------------------------------
            
            
        elseif k>N1 && k<N2
            
            
           %---------------------------Prediction from k-1 to k
            
            %............................Dynamic model coefficients
            %............. CN|k
            
            CNk = zeros(4,4);
            
            for ii=0:N2-k-1
                
                CNk = CNk + F^ii*Q*(F^ii)';
                
            end
            %.............
            
            Gk = Q - Q*F^(N2-k)'/(CNk + F^(N2-k)*Q*F^(N2-k)')*F^(N2-k)*Q;
            
            Gk_N = Gk*F^(N2-k)'/(CNk);
            Gk_km1 = F - Gk_N*F^(N2-k+1);
            
            
            
            Gy1 = horzcat(Gk_km1,Gk_N);
            Gy2 = horzcat(zeros(4,4),eye(4));
            Gy = vertcat(Gy1,Gy2);
            
            GCy = blkdiag(Gk,zeros(4,4));
            
            %.............
            
            Yp = Gy*Ys;
            Pyp = Gy*Pys*Gy' + GCy;
            
            %--------------------------------------------------
            
            Pyh = Pyp;
            
            Yh(:,k+1) = Yp ;
            
            Xh(:,k+1) = [eye(4) zeros(4,4)]*Yh(:,k+1);
            Ph = [eye(4) zeros(4,4)]*Pyh*[eye(4) zeros(4,4)]';
            
            %..................
            %.................. For the next time
            
            Ys = Yh(:,k+1);
            Pys = Pyh;
            
                
            %+++++++++++++++++++++
            
            %--------------------------------------------------------
            
        elseif k == N2

%             %+++++++++++++++++++++
%             %------------------------- N3
%             G_N3_N2 = C_N3_N2/XN2_Cov_t;
%             Xh(:,N3+1) = XN3_m_t + G_N3_N2*(Xh(:,N2+1) - XN2_m_t);
%             GN3 = XN3_Cov_t - C_N3_N2/XN2_Cov_t*C_N3_N2';
%             
%             Ph_N3 = GN3 + G_N3_N2*Ph_N2*G_N3_N2';
%             
%             
%             %-------------------N2 N3
%             
%             P_N3_N2 = G_N3_N2*Ph_N2;
%             %-------------------------- y 
%             
%             Ys = [Xh(:,N2+1);Xh(:,N3+1)];
%             Pys = [Ph_N2 P_N3_N2';
%                    P_N3_N2 Ph_N3];
            
               Xs = Xh(:,N2+1);
               Ps = Ph_N2;
               
               %-------------------------------
            %.........................
            

        elseif k>N2
            
          %---------------------------Prediction from k-1 to k
            
%             %............................Dynamic model coefficients
%             %............. CN|k
%             
%             CNk = zeros(4,4);
%             
%             for ii=0:N3-k-1
%                 
%                 CNk = CNk + F^ii*Q*(F^ii)';
%                 
%             end
%             %.............
%             
%             Gk = Q - Q*F^(N3-k)'/(CNk + F^(N3-k)*Q*F^(N3-k)')*F^(N3-k)*Q;
%             
%             Gk_N = Gk*F^(N3-k)'/(CNk);
%             Gk_km1 = F - Gk_N*F^(N3-k+1);
%             
%             
%             
%             Gy1 = horzcat(Gk_km1,Gk_N);
%             Gy2 = horzcat(zeros(4,4),eye(4));
%             Gy = vertcat(Gy1,Gy2);
%             
%             GCy = blkdiag(Gk,zeros(4,4));
%             
%             %.............
%             
%             Yp = Gy*Ys;
%             Pyp = Gy*Pys*Gy' + GCy;
%             
%             %--------------------------------------------------
%             
%             Pyh = Pyp;
%             
%             Yh(:,k+1) = Yp ;
%             
%             Xh(:,k+1) = [eye(4) zeros(4,4)]*Yh(:,k+1);
%             Ph = [eye(4) zeros(4,4)]*Pyh*[eye(4) zeros(4,4)]';
%             
%             %..................
%             %.................. For the next time
%             
%             Ys = Yh(:,k+1);
%             Pys = Pyh;
%             
%             
%             
%             %+++++++++++++++++++++

Xp = F*Xs;
Pp = F*Ps*F' + Q;

Xh(:,k+1) = Xp;
Ph = Pp;

Xs = Xp;
Ps = Pp;

            %--------------------------------------------------------
            
            
        end

    end

    %-----------------------


    %-----------------------

    EEp(iter,:) = sqrt((Xr(1,:) - Xh(1,:)).^2 + (Xr(3,:) - Xh(3,:)).^2);
    EEv(iter,:) = sqrt((Xr(2,:) - Xh(2,:)).^2 + (Xr(4,:) - Xh(4,:)).^2);
%     EEz(iter,:) = sqrt((Z(1,:) - Xr(1,:)).^2 + (Z(2,:) - Xr(3,:)).^2);


    %==========================================================

end  % for iter=1:iteration

TT = toc/iteration
%-----------------------------

AEEp = mean(EEp,1);
AEEv = mean(EEv,1);
% AEEz = mean(EEz,1);

%----------------------------

kk = 0:N3;

figure(1)
hold on
plot(kk,log10(AEEp),'o')
% hold on
% plot(kk,AEEz,'k')

% figure(2)
% hold on
% plot(kk,AEEv,'k')


%------------------ Test

C = [XN0_Cov_r C_N1_N0';C_N1_N0 XN0_Cov_r]

eig(C)

%-------------------
















