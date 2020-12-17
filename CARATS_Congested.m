clear;
close all;
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

A1 = signals(1).data;
A2 = signals(2).data;
A3 = signals(3).data;

%%1st table
data_FLT_1                 = table2array(A1(:,2));
%data_FLT_num        = table2array(A1(:,7));
data_time_ALL_1       = table2array(A1(:,1));
%[hh,mm,ss]              = hms(data_time);
%data_time_sec         = hh*60^2 + mm*60 + ss;
data_LAT_ALL_1        = table2array(A1(:,3));
data_LON_ALL_1       = table2array(A1(:,4));
data_ALT_ALL_1        = table2array(A1(:,5));
data_FLT_type_1        = table2array(A1(:,6));

%%2st table
data_FLT_2                 = table2array(A2(:,2));
%data_FLT_num        = table2array(A2(:,7));
data_time_ALL_2       = table2array(A2(:,1));
%[hh,mm,ss]              = hms(data_time);
%data_time_sec         = hh*60^2 + mm*60 + ss;
data_LAT_ALL_2        = table2array(A2(:,3));
data_LON_ALL_2       = table2array(A2(:,4));
data_ALT_ALL_2        = table2array(A2(:,5));
data_FLT_type_2        = table2array(A2(:,6));

%%3st table
data_FLT_3                 = table2array(A3(:,2));
%data_FLT_num        = table2array(A3(:,7));
data_time_ALL_3       = table2array(A3(:,1));
%[hh,mm,ss]              = hms(data_time);
%data_time_sec         = hh*60^2 + mm*60 + ss;
data_LAT_ALL_3        = table2array(A3(:,3));
data_LON_ALL_3       = table2array(A3(:,4));
data_ALT_ALL_3        = table2array(A3(:,5));
data_FLT_type_3        = table2array(A3(:,6));

%%After Vertcat of all
data_FLT_ALL                 = [data_FLT_1;data_FLT_2;data_FLT_3];
%data_FLT_num        = table2array(A1(:,7));
data_time_ALL       = [data_time_ALL_1; data_time_ALL_2; data_time_ALL_3];
%[hh,mm,ss]              = hms(data_time);
%data_time_sec         = hh*60^2 + mm*60 + ss;
data_LAT_ALL        = [data_LAT_ALL_1; data_LAT_ALL_2; data_LAT_ALL_3];
data_LON_ALL       = [data_LON_ALL_1; data_LON_ALL_2; data_LON_ALL_3];
data_ALT_ALL        = [data_ALT_ALL_1; data_ALT_ALL_2; data_ALT_ALL_3];
data_FLT_type        = [data_FLT_type_1; data_FLT_type_2; data_FLT_type_3];

%%Serch by Flight names
% k = contains(data_FLT,"AP01162");
% %data_time  = data_time_sec(k);
% data_FLTtype = data_FLT_type(k);
% data_LAT   = data_LAT_ALL(k);
% data_LON   = data_LON_ALL(k);
% data_ALT   = data_ALT_ALL(k);

%%Search by Flight distenations
k = find((data_ALT_ALL<50) & (abs(139.7 - data_LON_ALL) < 0.1) & (abs(35.5 - data_LAT_ALL)<0.1) );
%data_time  = data_time_sec(k);
data_FLT_PICKED = data_FLT_ALL(k);
 f = contains(data_FLT_ALL,data_FLT_PICKED);
%data_time  = data_time_sec(f);
data_FLT = data_FLT_ALL(f);
data_FLTtype = data_FLT_type(f);
data_LAT   = data_LAT_ALL(f);
data_LON   = data_LON_ALL(f);
data_ALT   = data_ALT_ALL(f);

figure(1);
for ii = 1:1:length(data_FLT)-1
    if isequal(data_FLT{ii},  data_FLT{ii+1})
        geoplot([data_LAT(ii),data_LAT(ii+1)],[data_LON(ii),data_LON(ii+1)],'-r');
        hold on;
        geolimits([34,36],[136.5,140.5]);
    end
end