clear;
close all;
rad = pi/180;
cd('C:\Users\Harada\Documents\MATLAB\data_library\CARATS_Open_Data_2017\20170417'); %Change here to observe the flights depending on the date
dirName = pwd;
files = dir('*.csv');
m = 1;
for n = 1:numel(files)
    if strfind(files(n).name, '.csv') > 0
        sprintf('%s', files(n).name)
        signals(m).filename = files(n).name;
        signals(m).data = readtable([dirName '\' files(n).name]);
        m = m+1;
    end
end
cd('../..')
cd('../data_library/map');
load('japanmap')
cd('../..')

A1 = signals(1).data;
A2 = signals(2).data;
A3 = signals(3).data;
DATA_table = [A1;A2;A3];
%% 
% data3D = table(zeros());
%  for m = 1:1:4276
%    for k =1:1:height(DATA(:,1))
%     if isequal(DATA(k,2), data_FLT_DATA(m,1))==1
%         data3D(k,1,m) = DATA(k,1);
%         data3D(k,2,m) = DATA(k,2);
%         data3D(k,3,m) = DATA(k,3);
%         data3D(k,4,m) = DATA(k,4);
%         data3D(k,5,m) = DATA(k,5);
%         data3D(k,6,m) = DATA(k,6);
%     end
%     end
%     
%end
%%1st table
data_FLT_1                 = table2array(A1(:,2));
%data_FLT_num        = table2array(A1(:,7));
data_time_ALL_1       = table2array(A1(:,1));
%[hh,mm,ss]                  = hms(data_time_ALL_1);
%data_time_sec_1         = hh*60^2 + mm*60 + ss;
data_LAT_ALL_1        = table2array(A1(:,3));
data_LON_ALL_1       = table2array(A1(:,4));
data_ALT_ALL_1        = table2array(A1(:,5));
data_FLT_type_1        = table2array(A1(:,6));

%%2nd table
data_FLT_2                 = table2array(A2(:,2));
%data_FLT_num        = table2array(A2(:,7));
data_time_ALL_2       = table2array(A2(:,1));
%[hh,mm,ss]                 = hms(data_time_ALL_2);
%data_time_sec_2         = hh*60^2 + mm*60 + ss;
data_LAT_ALL_2        = table2array(A2(:,3));
data_LON_ALL_2       = table2array(A2(:,4));
data_ALT_ALL_2        = table2array(A2(:,5));
data_FLT_type_2        = table2array(A2(:,6));

%%3rd table
data_FLT_3                 = table2array(A3(:,2));
%data_FLT_num        = table2array(A3(:,7));
data_time_ALL_3       = table2array(A3(:,1));
%[hh,mm,ss]                  = hms(data_time_ALL_3);
%data_time_sec_3         = hh*60^2 + mm*60 + ss;
data_LAT_ALL_3        = table2array(A3(:,3));
data_LON_ALL_3       = table2array(A3(:,4));
data_ALT_ALL_3        = table2array(A3(:,5));
data_FLT_type_3        = table2array(A3(:,6));

%%After Vertcat of all
data_FLT_ALL                 = [data_FLT_1;data_FLT_2;data_FLT_3];
%data_FLT_num        = cell2array(A1(:,7));
data_time_ALL       = [data_time_ALL_1; data_time_ALL_2; data_time_ALL_3];
%[hh,mm,ss]              = hms(data_time);
%data_time_sec         = hh*60^2 + mm*60 + ss;
data_LAT_ALL        = [data_LAT_ALL_1; data_LAT_ALL_2; data_LAT_ALL_3];
data_LON_ALL       = [data_LON_ALL_1; data_LON_ALL_2; data_LON_ALL_3];
data_ALT_ALL        = [data_ALT_ALL_1; data_ALT_ALL_2; data_ALT_ALL_3];
data_FLT_type        = [data_FLT_type_1; data_FLT_type_2; data_FLT_type_3];
DATA = [data_time_ALL,data_FLT_ALL,data_LAT_ALL,data_LON_ALL,data_ALT_ALL,data_FLT_type];


%DATA = [data_time_ALL,data_FLT_ALL,data_LAT_ALL,data_LON_ALL,data_ALT_ALL,data_FLT_type];
%%Serch by Flight names
% k = contains(data_FLT,"AP01162");
% %data_time  = data_time_sec(k);
% data_FLTtype = data_FLT_type(k);
% data_LAT   = data_LAT_ALL(k);
% data_LON   = data_LON_ALL(k);
% data_ALT   = data_ALT_ALL(k);
%%Search by Flight distenations
k = find((data_ALT_ALL<50) & (abs(data_LON_ALL - 139.7) < 0.1) & (abs(data_LAT_ALL - 35.5)<0.1));
%data_time  = data_time_sec(k);
data_FLT_PICKED = data_FLT_ALL(k);
f = contains(data_FLT_ALL,data_FLT_PICKED);
%data_time  = data_time_sec(f);
data_FLT = data_FLT_ALL(f);
data_FLTtype = data_FLT_type(f);
data_LAT   = data_LAT_ALL(f);
data_LON   = data_LON_ALL(f);
data_ALT   = data_ALT_ALL(f);

data_FLT_DATA = unique(data_FLT_PICKED);
sample = contains(data_FLT_ALL,data_FLT_DATA(1));
first = DATA_table(sample);
% for m = 1:1:300
%    for k =1:1:length(data_FLT(:,1))
%     if isequal(data_FLT(k), data_FLT_DATA(m,1))==1
%         data3D(k,:,m) = data_FLT(k,:);
%     end
%    end
% end
%d = find((data_ALT_ALL<1000) & (abs(data_LON_ALL - 139.7) < 0.1) & (abs(data_LAT_ALL - 35.5)<0.1));
%Calculate distances between airplane positions and within150[NM] from
%ARLON
ARLON = [139 + 58/60 + 59.8/3600, 35 + 15/60 + 25.3/3600]; %degree
TIA = [139 + 46/60 + 52/3600, 35 + 33/60 + 12/3600]; %degree
FUK = [130 + 27/60 + 6/3600 , 33 + 35/60 +4/3600];
a = 6378137;
b = 6356752.314;
e2 = 1 - (b^2/a^2);
L = length(data_LON);
mu_y = zeros(L,1);
W = zeros(L,1);
M = zeros(L,1);
N = zeros(L,1);
dx = zeros(L,1);
dy = zeros(L,1);
D = zeros(L,1);
D_km = zeros(L,1);
% for i = 1:L
%     mu_y(i,1) = deg2rad((ARLON(2)) + deg2rad(data_LAT(i))/2);
%     W(i,1) = sqrt(1- e2 * ((sin(mu_y(i,1)).^2)));
%     M(i,1) = a*(1 - e2)/(W(i,1).^3);
%     N(i,1) = a/W(i,1);
%     dy(i,1) = deg2rad(data_LAT(i)) - deg2rad(ARLON(2));
%     dx(i,1) = deg2rad(data_LON(i)) - deg2rad(ARLON(1));
%     D(i,1) = sqrt((dy(i,1)*M(i,1)).^2 + (dx(i,1)*N(i,1)*cos(mu_y(i,1))).^2);
% end

% for i = 1:L
%     mu_y(i,1) = deg2rad((TIA(2)) + deg2rad(data_LAT(i))/2);
%     W(i,1) = sqrt(1- e2 * ((sin(mu_y(i,1)).^2)));
%     M(i,1) = a*(1 - e2)/(W(i,1).^3);
%     N(i,1) = a/W(i,1);
%     dy(i,1) = deg2rad(data_LAT(i)) - deg2rad(TIA(2));
%     dx(i,1) = deg2rad(data_LON(i)) - deg2rad(TIA(1));
%     D(i,1) = sqrt((dy(i,1)*M(i,1)).^2 + (dx(i,1)*N(i,1)*cos(mu_y(i,1))).^2)/1000; %[km]
% end
for i = 1:L
    mu_y(i,1) = deg2rad((TIA(2)) + deg2rad(data_LAT(i))/2);
    W(i,1) = sqrt(1- e2 * ((sin(mu_y(i,1)).^2)));
    M(i,1) = a*(1 - e2)/(W(i,1).^3);
    N(i,1) = a/W(i,1);
    dy(i,1) = deg2rad(-data_LAT(i)) + deg2rad(TIA(2));
    dx(i,1) = deg2rad(-data_LON(i)) + deg2rad(TIA(1));
    D(i,1) = sqrt((dy(i,1)*M(i,1)).^2 + (dx(i,1)*N(i,1)*cos(mu_y(i,1))).^2)/1000; %[km]
end
diffALT = diff(data_ALT);
p = find(((diff(D)>-10)&(0>diff(D))));
Mem_FLT = unique(data_FLT_PICKED);
Lia_FLT = ismember(data_FLT_ALL,Mem_FLT);
data_FLT_FUK2TIA = data_FLT(p);
data_FLTtype_FUK2TIA = data_FLT_type(p);
data_LAT_FUK2TIA   = data_LAT(p);
data_LON_FUK2TIA   = data_LON(p);
data_ALT_FUK2TIA   = data_ALT(p);
D_FUK2TIA = D(p); 


% q = contains(data_FLT,data_FLT_FUK2TIA);
% data_FLT_FUK2TIA_final = data_FLT_ALL(q);
% g = zeros(L,1);
% for i = 3:L
%     g(i-2) = find((D(i)-D(i-2)) < 0);
% end
% data_FLT_FUK2TIA = data_FLT_ALL(g);
% p = contains(data_FLT_ALL,data_FLT_FUK2TIA);
% D_FUK2TIA = D(p);
%  for ii = 1:1:length(data_FLT_PICKED)
%     d(ii) = find(abs(data_LON(ii)-data_LON(ii+1))>10);
%     lon = data_LON(d);
%   
% end
% for m = 1:1:length(data_FLT(f))
%     lat(,m) = 
%     if data_FLT(m) = data_FLT(m+1)
%         data_FLT_();
%     end
%     if abs(data_LAT(m) - data_LAT(m+1))>1
%         lat(n+1,m) = lat(n,m); 
%     end
% end
% for m = 1:1:length(data_FLT)
%     for n=1:1:length(data_FLT)
%      if isequal(data_FLT{m},  data_FLT{m+1})
%          lat(m,n) = data_LAT(m);
%          lon(m,n) = data_LON(m); 
%          alt(m,n) = data_ALT(m);
%      else
%          break;
%      end
%     end
% end
% figure(1);

% 
% for ii = 1:1:length(data_FLT)-1
%    
%         geoplot([data_LAT(ii),data_LAT(ii+1)],[data_LON(ii),data_LON(ii+1)],'-r');
%         hold on;
%         geolimits([34,36],[136.5,140.5]);
%     end
% end
