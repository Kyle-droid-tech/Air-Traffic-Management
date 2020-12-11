clear
close all
cdir = cd;

rad = pi/180;
deg = 180/pi;
ml2km = 1.852;
ft2m  = 0.3048;
kt2ms = 0.5144;
R0 = 6378.137*10^3; %Mean radius of the Earth [m]
g0 = 9.80665;

cd('../data_library/map');
load('japanmap')
cd('../..')

%% Load Aircraft Performance Model
ACtype = 'B772';
cd('data_library/BADA_3_14');
OPF = load(['OPF_',ACtype,'_BADA3-14.data']);
cd(cdir)

%% Concentric circle
ARLON = [139+58/60+59.8/3600, 35+15/60+25.3/3600]*rad;

dlon_100NM = 100*ml2km*10^3/R0/cos(ARLON(1,2))*deg;
dlat_100NM = 100*ml2km*10^3/R0*deg;

theta = 0 : 2*pi/720 : 2*pi;

x_50NM  = ARLON(1,1)*deg+0.5*dlon_100NM*cos(theta);
y_50NM  = ARLON(1,2)*deg+0.5*dlat_100NM*sin(theta);
x_100NM = ARLON(1,1)*deg+1.0*dlon_100NM*cos(theta);
y_100NM = ARLON(1,2)*deg+1.0*dlat_100NM*sin(theta);
x_150NM = ARLON(1,1)*deg+1.5*dlon_100NM*cos(theta);
y_150NM = ARLON(1,2)*deg+1.5*dlat_100NM*sin(theta);

%%
%Initial position
pos_A(1,:) = [137.102473463831, 34.4240061313947]*rad;
pos_B(1,:) = [137.803514855962, 33.5078947891613]*rad;

r0_A = [cos(pos_A(1,2))*cos(pos_A(1,1)); cos(pos_A(1,2))*sin(pos_A(1,1)); sin(pos_A(1,2))];
r0_B = [cos(pos_B(1,2))*cos(pos_B(1,1)); cos(pos_B(1,2))*sin(pos_B(1,1)); sin(pos_B(1,2))];
rf = [cos(ARLON(1,2))*cos(ARLON(1,1)); cos(ARLON(1,2))*sin(ARLON(1,1)); sin(ARLON(1,2))];

r2_A = cross(r0_A,rf)/norm(cross(r0_A,rf),2);
r1_A = cross(r2_A,r0_A)/norm(cross(r2_A,r0_A),2);
r2_B = cross(r0_B,rf)/norm(cross(r0_B,rf),2);
r1_B = cross(r2_B,r0_B)/norm(cross(r2_B,r0_B),2);

x = (0 : 10 : 150) *ml2km*10^3;

for ii = 2 : size(x,2)
   [pos_A(ii,2), pos_A(ii,1)] = eta_zeta2phi_theta(x(ii)/R0,0,r0_A,r1_A,r2_A);
   [pos_B(ii,2), pos_B(ii,1)] = eta_zeta2phi_theta(x(ii)/R0,0,r0_B,r1_B,r2_B);
end

alt = 35000 : (5200-35000)/(size(x,2)-1) : 5200;

Cons_pos_A = pos_A * deg;
Cons_pos_B = pos_B * deg;

Cons_alt = linspace(35000 , 5200, (35000-5200));
%% CAS optimization by DP
CAS_A = 200 : 10 : 300;
CAS_B = 200 : 10 : 300;

a_A = 0;   a_B = 0;   % Case 1
% a_A = 0.5; a_B = 0.5; % Case 2
% a_A = 0;   a_B = 0.3; % Case 3
 %a_A = 0;   a_B = 0.3; % Case 4
t_sep = 90;

initialization

CAS_0 = 280; % Case 1 to 3
CAS_f = 210;
alt_f = alt(size(x,2));

kk = size(x,2)-1;
N = kk;
[T_f, rho_f, p_f, a_f] = ISA(alt_f*ft2m);
TAS_f = cas2tas(rho_f,p_f*100,CAS_f*kt2ms);

[~,rho_A,p_A,~] = ISA(alt(kk)*ft2m);
[~,rho_B,p_B,~] = ISA(alt(kk)*ft2m);

for ii = 1 : size(CAS_A,2) % CAS of A
    TAS_A_1 = cas2tas(rho_A,p_A*100,CAS_A(ii)*kt2ms);
    [dt_A, FF_A] = FF_BADA(x(kk),x(kk+1),alt(kk),alt(kk+1),TAS_A_1,TAS_f,OPF);
    dF_A = FF_A*dt_A;
    
    for jj = 1 : size(CAS_B,2) % CAS of B
        TAS_B_1 = cas2tas(rho_B,p_B*100,CAS_B(jj)*kt2ms);
        [dt_B, FF_B] = FF_BADA(x(kk),x(kk+1),alt(kk),alt(kk+1),TAS_B_1,TAS_f,OPF);
        dF_B = FF_B*dt_B;
        
        J_opt(ii,jj,kk) = dF_A+a_A*dt_A + dF_B+a_B*dt_B;
        F_opt(ii,jj,kk) = dF_A+dF_B;
        F_opt_A(ii,jj,kk) = dF_A;
        F_opt_B(ii,jj,kk) = dF_B;
        T_opt_A(ii,jj,kk) = dt_A;
        T_opt_B(ii,jj,kk) = dt_B;
        ind_CAS_A(ii,jj,kk+1) = 2;
        ind_CAS_B(ii,jj,kk+1) = 2;
    end
end

%% Optimization from kk=1 to kk=N
for kk = size(x,2)-2 : -1 : 1
    % disp(kk)
    [~,rho_A_1,p_A_1,~] = ISA(alt(kk)  *ft2m);
    [~,rho_A_2,p_A_2,~] = ISA(alt(kk+1)*ft2m);
    [~,rho_B_1,p_B_1,~] = ISA(alt(kk)  *ft2m);
    [~,rho_B_2,p_B_2,~] = ISA(alt(kk+1)*ft2m);
    
    for ii = 1 : size(CAS_A,2) % CAS of A at PREVIOUS point
        TAS_A_1 = cas2tas(rho_A_1,p_A_1*100,CAS_A(ii)*kt2ms);
        
        for jj = 1 : size(CAS_B,2) % CAS of B at PREVIOUS point
            TAS_B_1 = cas2tas(rho_B_1,p_B_1*100,CAS_B(jj)*kt2ms);
            
            for mm = 1 : size(CAS_A,2) % CAS of A at NEXT point
                TAS_A_2 = cas2tas(rho_A_2,p_A_2*100,CAS_A(mm)*kt2ms);
                [dt_A, FF_A] = FF_BADA(x(kk),x(kk+1),alt(kk),alt(kk+1),TAS_A_1,TAS_A_2,OPF);
                dF_A = FF_A*dt_A;

                for nn = 1 : size(CAS_B,2) % CAS of B at NEXT point
                    TAS_B_2 = cas2tas(rho_B_2,p_B_2*100,CAS_B(nn)*kt2ms);
                    [dt_B, FF_B] = FF_BADA(x(kk),x(kk+1),alt(kk),alt(kk+1),TAS_B_1,TAS_B_2,OPF);
                    dF_B = FF_B*dt_B;

                    J_all(mm,nn) = J_opt(mm,nn,kk+1)+ dF_A+a_A*dt_A + dF_B+a_B*dt_B;
                    F_all(mm,nn) = F_opt(mm,nn,kk+1)+dF_A+dF_B;
                    F_all_A(mm,nn) = F_opt_A(mm,nn,kk+1)+dF_A;
                    F_all_B(mm,nn) = F_opt_B(mm,nn,kk+1)+dF_B;
                    T_all_A(mm,nn) = T_opt_A(mm,nn,kk+1)+dt_A;
                    T_all_B(mm,nn) = T_opt_B(mm,nn,kk+1)+dt_B;
                end
            end
            J_opt(ii,jj,kk) = min(min(J_all));
            [minval1, min_ind1]   = min(J_all); 
            [minval2, min_ind2]   = min(minval1);
            ind_CAS_B(ii,jj,kk+1) = min_ind2;
            ind_CAS_A(ii,jj,kk+1) = min_ind1(min_ind2);
            F_opt(ii,jj,kk) = F_all(ind_CAS_A(ii,jj,kk+1),ind_CAS_B(ii,jj,kk+1));
            F_opt_A(ii,jj,kk) = F_all_A(ind_CAS_A(ii,jj,kk+1),ind_CAS_B(ii,jj,kk+1));
            F_opt_B(ii,jj,kk) = F_all_B(ind_CAS_A(ii,jj,kk+1),ind_CAS_B(ii,jj,kk+1));
            T_opt_A(ii,jj,kk) = T_all_A(ind_CAS_A(ii,jj,kk+1),ind_CAS_B(ii,jj,kk+1));
            T_opt_B(ii,jj,kk) = T_all_B(ind_CAS_A(ii,jj,kk+1),ind_CAS_B(ii,jj,kk+1));
        end
    end
end

% J_opt_0 = J_opt(:,:,1);                            % Case 4
 %delt_f = T_opt_A(:,:,1)-T_opt_B(:,:,1);            % Case 4
 %J_opt_0 = J_opt_0.*(delt_f>=90)+10^6.*(delt_f<90); % Case 4


%%
ind_CAS_A_opt(1) = find(abs(CAS_A-CAS_0)<10^-6); % Case 1 to 3 
ind_CAS_B_opt(1) = find(abs(CAS_B-CAS_0)<10^-6); % Case 1 to 3

% [val, ind] = min(J_opt_0,[],1);           % Case 4
% [~, ind_CAS_B_opt(1)] = min(val);         % Case 4
% ind_CAS_A_opt(1) = ind(ind_CAS_B_opt(1)); % Case 4


for kk = 2 : size(x,2)
    ind_CAS_A_opt(kk) = ind_CAS_A(ind_CAS_A_opt(kk-1),ind_CAS_B_opt(kk-1),kk);
    ind_CAS_B_opt(kk) = ind_CAS_B(ind_CAS_A_opt(kk-1),ind_CAS_B_opt(kk-1),kk);
end

CAS_A_opt=CAS_A(ind_CAS_A_opt);
CAS_B_opt=CAS_B(ind_CAS_B_opt);
 
for kk = 1 : size(x,2)
    [~,rho_A(kk,1),p_A(kk,1),~] = ISA(alt(kk)*ft2m);
    TAS_A_opt(kk,1) = cas2tas(rho_A(kk,1),p_A(kk,1)*100,CAS_A_opt(kk)*kt2ms);
%     [~,rho_B(kk,1),p_B(kk,1),~] = ISA(alt(kk)*ft2m);
    rho_B(kk,1) = rho_A(kk,1);
    p_B(kk,1)   = p_A(kk,1);
    TAS_B_opt(kk,1) = cas2tas(rho_B(kk,1),p_B(kk,1)*100,CAS_B_opt(kk)*kt2ms);
end


dt_opt_A = diff(x')./((TAS_A_opt(1:size(x,2)-1)+TAS_A_opt(2:size(x,2)))/2);
dt_opt_B = diff(x')./((TAS_B_opt(1:size(x,2)-1)+TAS_B_opt(2:size(x,2)))/2);
%    dt_opt_A = diff(x(kk)')./((TAS_A_opt(kk)+TAS_A_opt(kk))/2);
%    dt_opt_B = diff(x(kk)')./((TAS_B_opt(kk)+TAS_B_opt(kk))/2);

time_A_opt(1,1) = 0;
time_B_opt(1,1) = 0;
for kk = 1 : size(x,2)-1
    time_A_opt(kk+1,1) = time_A_opt(kk,1) + dt_opt_A(kk,1);
    time_B_opt(kk+1,1) = time_B_opt(kk,1) + dt_opt_B(kk,1);
end

%%Protection area
%Calcurate distance between k and k+1 till k=16 for airplane A
a_A = zeros(size(x,2)-1,1);
for ii = 1 : 1 :size(x,2)-1
    a_A(ii) = (TAS_A_opt(ii+1) - TAS_A_opt(ii)) / dt_opt_A(ii);
end
x_det_A = zeros(2000);
for jj = 1 : 1  : size(dt_opt_A,1)
     for ii = 1 : 1 : size(TAS_A_opt,1)
        for t = 1 : 1 : fix(dt_opt_A(jj))
              x_det_A(t,jj) = x(jj)+  + TAS_A_opt(jj)*t + (a_A(jj)*t^2)/2;
        end
     end     
end
  x_det_A(x_det_A==0) = [];
 x_det_A =[0,(x_det_A)];
%Calcurate distance between k and k+1 till k=16 for airplane B 
a_B = zeros(size(x,2)-1,1);
for ii = 1 : 1 :size(x,2)-1
    a_B(ii) = (TAS_B_opt(ii) - TAS_B_opt(ii+1)) / dt_opt_B(ii);
end
x_det_B = zeros(2000);
for jj = 1 : 1  : size(dt_opt_B,1)
    for ii = 1 : 1 : size(TAS_B_opt,1)
        for t = 1 : 1 : fix(dt_opt_B(jj))
              x_det_B(t,jj) =x(jj) +  TAS_B_opt(jj)*t + (a_B(jj)*t^2)/2;
         end
    end     
end
  x_det_B(x_det_B==0) = [];
  x_det_B =[0, (x_det_B)];
  
  %
  pos_det_A(1,:) = pos_A(1,:);
  pos_det_B(1,:) = pos_B(1,:);
for ii = 2 : size(x_det_A,2)
   [pos_det_A(ii,2), pos_det_A(ii,1)] = eta_zeta2phi_theta(x_det_A(ii)/R0,0,r0_A,r1_A,r2_A);
   [pos_det_B(ii,2), pos_det_B(ii,1)] = eta_zeta2phi_theta(x_det_B(ii)/R0,0,r0_B,r1_B,r2_B);
end
p_det_A = pos_det_A / rad;
p_det_B = pos_det_B / rad;

%%Distance Calcuration between airplane A and B every 1 sec
a = 6378137;
b = 6356752.314;
e2 = 1 - (b^2/a^2);

L = length(x_det_A)-1;

mu_y = zeros(L,1);
W = zeros(L,1);
M = zeros(L,1);
N = zeros(L,1);
dx = zeros(L,1);
dy = zeros(L,1);
D = zeros(L,1);
%for A
for i=1:L
    mu_y(i,1) = deg2rad((pos_det_A(i,2) + pos_det_B(i,2))/2);
    W(i,1) = sqrt(1-e2 * ((sin(mu_y(i,1)).^2)));
    M(i,1) = a*(1 - e2)/(W(i,1).^3);
    N(i,1) = a/W(i,1);
    dy(i,1) = deg2rad(pos_det_B(i,2) - pos_det_A(i,2));
    dx(i,1) = deg2rad(pos_det_B(i,1) - pos_det_A(i,1));
    
    D(i,1) = sqrt((dy(i,1)*M(i,1)).^2 + (dx(i,1)*N(i,1)*cos(mu_y(i,1))).^2);
end

%%
% disp(J_opt(ind_CAS_A_opt(1),ind_CAS_B_opt(1),1));
disp(['Topt(A) [s] = ',num2str(T_opt_A(ind_CAS_A_opt(1),ind_CAS_B_opt(1),1))]);
disp(['Topt(B) [s] = ',num2str(T_opt_B(ind_CAS_A_opt(1),ind_CAS_B_opt(1),1))]);
disp(['Fopt(A) [s] = ',num2str(F_opt_A(ind_CAS_A_opt(1),ind_CAS_B_opt(1),1))]);
disp(['Fopt(B) [s] = ',num2str(F_opt_B(ind_CAS_A_opt(1),ind_CAS_B_opt(1),1))]);

%%
fs  = 14;
 figure;
 hold on; grid on;
 plot(lon,lat,'k')
 axis equal
 plot(x_50NM,y_50NM,'r-.'); plot(x_100NM,y_100NM,'r-.'); plot(x_150NM,y_150NM,'r-.');
 plot(pos_A(1,1)*deg,pos_A(1,2)*deg,'ro'); plot(pos_B(1,1)*deg,pos_B(1,2)*deg,'bo');
 plot([pos_A(1:15,1); ARLON(1,1) ]*deg, [pos_A(1:15,2); ARLON(1,2)]*deg,'-r.')
 plot([pos_B(1:15,1); ARLON(1,1) ]*deg, [pos_B(1:15,2); ARLON(1,2)]*deg,'-b.')
 xlabel('Longitude [deg]','fontweight','bold','fontsize',fs)
 ylabel('Latitude [deg]','fontweight','bold','fontsize',fs)
 set(gca,'fontweight','bold','fontsize',fs)

%%
x2 = -flip(x)/ml2km/10^3;

figure;
hold on; grid on;
p1 = plot(x2,CAS_A(ind_CAS_A_opt),'-r');
plot(x2(2:end-1),CAS_A(ind_CAS_A_opt(2:end-1)),'ro','MarkerFaceColor','r','MarkerSize',3)
p2 = plot(x2,CAS_B(ind_CAS_B_opt),'--b');
plot(x2(2:end-1),CAS_B(ind_CAS_B_opt(2:end-1)),'bo','MarkerFaceColor','b','MarkerSize',3)
plot(x2(1,1),CAS_A(ind_CAS_A_opt(1)),'ro')
plot(x2(1,1),CAS_B(ind_CAS_B_opt(1)),'bo')
p3 = plot(x2(1,end),210,'mo');
legend([p1,p2,p3],'Aircraft A','Aircraft B','Terminal point')
set(gca,'YLim',[200 300])
xlabel('$$ x \ [\mathrm{NM}] $$','Interpreter','latex','fontweight','bold','fontsize',fs)
ylabel('$$ V_{\mathrm{CAS}} \ [\mathrm{kt}] $$','Interpreter','latex','fontweight','bold','fontsize',fs)
 text(-140,290,'$$ J=\Sigma \int \mu \mathrm{d}t $$','Interpreter','latex','fontweight','bold','fontsize',fs-2)
% text(-140,290,'$$ J=\Sigma \int (\mu+0.5) \mathrm{d}t $$','Interpreter','latex','fontweight','bold','fontsize',fs-2)
%text(-143,290,'$$ J=\int \mu_A \mathrm{d}t_A+\int (\mu_B+0.3) \mathrm{d}t_B $$','Interpreter','latex','fontweight','bold','fontsize',fs-2)
set(gca,'fontweight','bold','fontsize',fs)

%%
 figure;
 hold on; grid on;
 plot(pos_A(:,1)*deg,CAS_A(ind_CAS_A_opt),'r')
 plot(pos_B(:,1)*deg,CAS_B(ind_CAS_B_opt),'b')
 plot(pos_A(1,1)*deg,CAS_A(ind_CAS_A_opt(1)),'ro')
 plot(pos_B(1,1)*deg,CAS_B(ind_CAS_B_opt(1)),'bo')
 plot(ARLON(1,1)*deg,210,'mo')
 legend('Aircraft A','Aircraft B')
 set(gca,'YLim',[200 300])
 title('$$ J=\Sigma \int \mu \mathrm{d}t $$','Interpreter','latex')
 title('$$ J=\Sigma \int (\mu+0.5) \mathrm{d}t $$','Interpreter','latex')