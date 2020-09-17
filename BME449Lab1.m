%%.......................................................................%%
%.........BME 449 Lab 1 - Biomechanics.................________...........%
%.........For Dr. Sarles...............................|/.||.\|...........%
%.........By Robert Borkoski..............................||..............%
%.........Last Revised September 6 2020.................__/\__............%
%%.......................................................................%%

%% initializers

clear all, close all, clc, format long g

%% import data

fprintf('importing data.'); %begin data aquisition

% bone dimensions (all samples)

data.dim = xlsread('Bone-sample-dimensions-Sect1.xlsx');

dot()
% Cold, fast samples
data.CF1 = xlsread('CF_Specimen_RawData_1.csv'); %specimen 1
dot();
data.CF2 = xlsread('CF_Specimen_RawData_2.csv'); %specimen 2
dot();
data.CF3 = xlsread('CF_Specimen_RawData_3.csv'); %specimen 3
dot();
data.CF4 = xlsread('CF_Specimen_RawData_4.csv'); %specimen 4
dot();

% Cold, medium samples (1)

data.CM1 = xlsread('CM_Specimen_RawData_1.csv'); %specimen 1
dot();
data.CM2 = xlsread('CM_Specimen_RawData_2.csv'); %specimen 2
dot();
data.CM3 = xlsread('CM_Specimen_RawData_3.csv'); %specimen 3
dot();
data.CM4 = xlsread('CM_Specimen_RawData_4.csv'); %specimen 4
dot();

% Cold, medium samples (2)

data.CMM1 = xlsread('CM2_Specimen_RawData_1.csv'); %specimen 1
dot();
data.CMM2 = xlsread('CM2_Specimen_RawData_2.csv'); %specimen 2
dot();
data.CMM3 = xlsread('CM2_Specimen_RawData_3.csv'); %specimen 3
dot();
data.CMM4 = xlsread('CM2_Specimen_RawData_4.csv'); %specimen 4
dot();

% room temperature, fast samples

data.RTF1 = xlsread('RTF_Specimen_RawData_1.csv'); %specimen 1
dot();
data.RTF2 = xlsread('RTF_Specimen_RawData_2.csv'); %specimen 2
dot();
data.RTF3 = xlsread('RTF_Specimen_RawData_3.csv'); %specimen 3
dot();
data.RTF4 = xlsread('RTF_Specimen_RawData_4.csv'); %specimen 4
dot();

% room temperature, medium samples

data.RTM1 = xlsread('RTM_Specimen_RawData_1.csv'); %specimen 1
dot();
data.RTM2 = xlsread('RTM_Specimen_RawData_2.csv'); %specimen 2
dot();
data.RTM3 = xlsread('RTM_Specimen_RawData_3.csv'); %specimen 3
dot();
data.RTM4 = xlsread('RTM_Specimen_RawData_4.csv'); %specimen 4
dot();

% room temperature, slow samples

data.RTS1 = xlsread('RTS_Specimen_RawData_1.csv'); %specimen 1
dot();
data.RTS2 = xlsread('RTS_Specimen_RawData_2.csv'); %specimen 2
dot();
data.RTS3 = xlsread('RTS_Specimen_RawData_3.csv'); %specimen 3
dot();
data.RTS4 = xlsread('RTS_Specimen_RawData_4.csv'); %specimen 4
dot();
data.RTS5 = xlsread('RTS_Specimen_RawData_5.csv'); %specimen 5
dot();

% warm, fast samples

data.WF1 = xlsread('WF_Specimen_RawData_1.csv'); %specimen 1
dot();
data.WF2 = xlsread('WF_Specimen_RawData_2.csv'); %specimen 2
dot();
data.WF3 = xlsread('WF_Specimen_RawData_3.csv'); %specimen 3
dot();
data.WF4 = xlsread('WF_Specimen_RawData_4.csv'); %specimen 4
dot();

% warm, medium samples

data.WM1 = xlsread('WM_Specimen_RawData_1.csv'); %specimen 1
dot();
data.WM2 = xlsread('WM_Specimen_RawData_2.csv'); %specimen 2
dot();
data.WM3 = xlsread('WM_Specimen_RawData_3.csv'); %specimen 3
dot();
data.WM4 = xlsread('WM_Specimen_RawData_4.csv'); %specimen 4
dot();

fprintf('complete'); % data aquisition complete

%% plotting - raw data

% displacement / force charts

figure(1)
subplot(2,2,1)
plot(data.WM1(:,2),data.WM1(:,3)); % plot force vs. displacement - cold fast spec 1
subplot(2,2,2)
plot(data.WM2(:,2),data.WM2(:,3)); % plot force vs. displacement - cold fast spec 2
subplot(2,2,3)
plot(data.WM3(:,2),data.WM3(:,3)); % plot force vs. displacement - cold fast spec 3
subplot(2,2,4)
plot(data.WM4(:,2),data.WM4(:,3)); % plot force vs. displacement - cold fast spec 4

%% Data processing

% obtain bone dimensions from raw experimental data

dim.dim_RTS = data.dim(1:5,:); % sample dimensions - room temp slow
dim.dim_CM = data.dim(11:14,:); % sample dimensions - cold med (1)
dim.dim_RTM = data.dim(20:23,:); % sample dimensions - room temp med
dim.dim_WM = data.dim(30:33,:); % sample dimensions - warm med
dim.dim_CMM = data.dim(39:42,:); % sample dimensions - cold med (2)
dim.dim_CF = data.dim(48:51,:); % sample dimensions - cold fast
dim.dim_RTF = data.dim(57:60,:); % sample dimensions - room temp fast
dim.dim_WF = data.dim(66:69,:); % sample dimensions - warm fast

% generate list of names of raw data elements (to streamline)

Data_names = {'CF1','CF2','CF3','CF4','CM1','CM2','CM3','CM4','CMM1','CMM2','CMM3','CMM4','RTF1','RTF2','RTF3','RTF4','RTM1','RTM2','RTM3','RTM4','RTS1','RTS2','RTS3','RTS4','RTS5','WF1','WF2','WF3','WF4','WM1','WM2','WM3','WM4'};

% % calculate slope of force / displacement curves
% for j = 1:length(Data_names)
%     Slope_names = strcat('slope_',Data_names{j});
%     for m = 1:length(data.(Data_names{j}))-1
%         Slope.(Slope_names)(m) = (data.(Data_names{j})(m+1,3) - data.(Data_names{j})(m,3))./(data.(Data_names{j})(m+1,2) - data.(Data_names{j})(m,2));
%     end
%     clear m
%     init_window = round(0.2.*length(data.(Data_names{j})));
%     init = mean(Slope.(Slope_names)(1:init_window)); % establish 10 pt moving average to compare slopes against
%     if init == -inf
%         for z = 1:length(data.(Data_names{j}))
%             if Slope.(Slope_names)(z) == -inf
%                 Slope.(Slope_names)(z) = 0;
%                 init = mean(Slope.(Slope_names)(1:init_window));
%                 if init ~= -inf
%                     break
%                 end
%                 clear z
%             end
%         end
%     end
%     lim = 0;
%     while lim == 0
%         for w = init_window:length(data.(Data_names{j}))
%             window = mean(Slope.(Slope_names)((w-round(init_window/10)):(w+round(init_window/10))));
%             if window < 0.5*init
%                 linreg(j) = w;
%                 clear w
%                 break
%             end
%         end
%     break    
%     end
% end
% clear init, clear init_window, clear j, clear lim, clear Slope_names, clear window

% calculate equivalent stiffness k

data.FD = xlsread('Inspection.xlsx'); %Force / displacement inspection data

for s = 1:length(Data_names)
    K.(Data_names{s}) = data.FD(s,5) ./ (data.FD(s4)./1000)
end
% for u = 1:length(Data_names)
%     K.(Data_names{u}) = (data.(Data_names{u})(linreg(u),3) - 0) / (data.(Data_names{u})(linreg(u),2) - 0);
% end
% clear u

% compute moment of inertia I
Dim_names = {'Dim_CF1','Dim_CF2','Dim_CF3','Dim_CF4','Dim_CM1','Dim_CM2','Dim_CM3','Dim_CM4','Dim_CMM1','Dim_CMM2','Dim_CMM3','Dim_CMM4','Dim_RTF1','Dim_RTF2','Dim_RTF3','Dim_RTF4','Dim_RTM1','Dim_RTM2','Dim_RTM3','Dim_RTM4','Dim_RTS1','Dim_RTS2','Dim_RTS3','Dim_RTS4','Dim_RTS5','Dim_WF1','Dim_WF2','Dim_WF3','Dim_WF4','Dim_WM1','Dim_WM2','Dim_WM3','Dim_WM4'};

dim.Dim_CF1 = dim.dim_CF(1,:)
dim.Dim_CF2 = dim.dim_CF(2,:)
dim.Dim_CF3 = dim.dim_CF(3,:)
dim.Dim_CF4 = dim.dim_CF(4,:)
dim.Dim_CM1 = dim.dim_CM(1,:)
dim.Dim_CM2 = dim.dim_CM(2,:)
dim.Dim_CM3 = dim.dim_CM(3,:)
dim.Dim_CM4 = dim.dim_CM(4,:)
dim.Dim_CMM1 = dim.dim_CMM(1,:)
dim.Dim_CMM2 = dim.dim_CMM(2,:)
dim.Dim_CMM3 = dim.dim_CMM(3,:)
dim.Dim_CMM4 = dim.dim_CMM(4,:)
dim.Dim_RTF1 = dim.dim_RTF(1,:)
dim.Dim_RTF2 = dim.dim_RTF(2,:)
dim.Dim_RTF3 = dim.dim_RTF(3,:)
dim.Dim_RTF4 = dim.dim_RTF(4,:)
dim.Dim_RTM1 = dim.dim_RTM(1,:)
dim.Dim_RTM2 = dim.dim_RTM(2,:)
dim.Dim_RTM3 = dim.dim_RTM(3,:)
dim.Dim_RTM4 = dim.dim_RTM(4,:)
dim.Dim_RTS1 = dim.dim_RTS(1,:)
dim.Dim_RTS2 = dim.dim_RTS(2,:)
dim.Dim_RTS3 = dim.dim_RTS(3,:)
dim.Dim_RTS4 = dim.dim_RTS(4,:)
dim.Dim_RTS5 = dim.dim_RTS(5,:)
dim.Dim_WF1 = dim.dim_WF(1,:)
dim.Dim_WF2 = dim.dim_WF(2,:)
dim.Dim_WF3 = dim.dim_WF(3,:)
dim.Dim_WF4 = dim.dim_WF(4,:)
dim.Dim_WM1 = dim.dim_WM(1,:)
dim.Dim_WM2 = dim.dim_WM(2,:)
dim.Dim_WM3 = dim.dim_WM(3,:)
dim.Dim_WM4 = dim.dim_WM(4,:)

% define radius
for r = 1:length(Dim_names)
    r0.(Data_names{r}) = (((dim.(Dim_names{r})(1,2) + dim.(Dim_names{r})(1,3))/2)/2)/1000
    ri.(Data_names{r}) = ((((dim.(Dim_names{r})(1,2) + dim.(Dim_names{r})(1,3))/2)/2) - dim.(Dim_names{r})(1,4))/1000
    I.(Data_names{r}) = (pi()./4)*((r0.(Data_names{r}).^4) - (ri.(Data_names{r}).^4))
end
clear q, clear r

% determine l

l = 2.5*2.54/100 % convert 2.5 in to m

% determine young's modulus

for o = 1:length(Data_names)
    E.(Data_names{o}) = (K.(Data_names{o})*(l.^3))/(48.*I.(Data_names{o}));
end
clear o

% calculate M

for y = 1:length(Data_names)
    M.(Data_names{y}) = ((data.(Data_names{y})(y,3))*l)./4
end
clear y,



% compute stresses

% for b = 1:length(Data_names)
%     for v = 1:length(M.(Data_names{b}))
%         Stress.(Data_names{b})(v) = (M.(Data_names{b})(v))*(r0.(Data_names{b}))/I.(Data_names{b});
%     end
% end
% clear b, clear v

% compute yield stress
for g = 1:length(Data_names)
    Yield_Points.(Data_names{g}) = (data.FD(g,5))./(pi().*(((r0.(Data_names{g})).^2 - ((ri.(Data_names{g}))).^2)));
end
clear g

% compute max stresses
for a = 1:length(Data_names)
    Max_Stress.(Data_names{a}) = max(data.(Data_names{a})(:,3))/(pi().*(r0.(Data_names{a}).^2 - ri.(Data_names{a}).^2));
end
clear a, clear z, clear var

%% prep for statistical testing

% form groups - cold fast

E_CF = [E.CF1,E.CF2,E.CF3,E.CF4];
Yield_Points_CF = [Yield_Points.CF1,Yield_Points.CF2,Yield_Points.CF3,Yield_Points.CF4];
Max_Stress_CF = [Max_Stress.CF1,Max_Stress.CF2,Max_Stress.CF3,Max_Stress.CF4];

% cold med (1)

E_CM = [E.CM1,E.CM2,E.CM3,E.CM4];
Yield_Points_CM = [Yield_Points.CM1,Yield_Points.CM2,Yield_Points.CM3,Yield_Points.CM4];
Max_Stress_CM = [Max_Stress.CM1,Max_Stress.CM2,Max_Stress.CM3,Max_Stress.CM4];

% cold med (2)

E_CMM = [E.CMM1,E.CMM2,E.CMM3,E.CMM4];
Yield_Points_CMM = [Yield_Points.CMM1,Yield_Points.CMM2,Yield_Points.CMM3,Yield_Points.CMM4];
Max_Stress_CMM = [Max_Stress.CMM1,Max_Stress.CMM2,Max_Stress.CMM3,Max_Stress.CMM4];

% room temp fast

E_RTF = [E.RTF1,E.RTF2,E.RTF3,E.RTF4];
Yield_Points_RTF = [Yield_Points.RTF1,Yield_Points.RTF2,Yield_Points.RTF3,Yield_Points.RTF4];
Max_Stress_RTF = [Max_Stress.RTF1,Max_Stress.RTF2,Max_Stress.RTF3,Max_Stress.RTF4];

% room temp med

E_RTM = [E.RTM1,E.RTM2,E.RTM3,E.RTM4];
Yield_Points_RTM = [Yield_Points.RTM1,Yield_Points.RTM2,Yield_Points.RTM3,Yield_Points.RTM4];
Max_Stress_RTM = [Max_Stress.RTM1,Max_Stress.RTM2,Max_Stress.RTM3,Max_Stress.RTM4];

% room temp slow

E_RTS = [E.RTS1,E.RTS2,E.RTS3,E.RTS4,E.RTS5];
Yield_Points_RTS = [Yield_Points.RTS1,Yield_Points.RTS2,Yield_Points.RTS3,Yield_Points.RTS4,Yield_Points.RTS5];
Max_Stress_RTS = [Max_Stress.RTS1,Max_Stress.RTS2,Max_Stress.RTS3,Max_Stress.RTS4,Max_Stress.RTS5];

% warm fast

E_WF = [E.WF1,E.WF2,E.WF3,E.WF4];
Yield_Points_WF = [Yield_Points.WF1,Yield_Points.WF2,Yield_Points.WF3,Yield_Points.WF4];
Max_Stress_WF = [Max_Stress.WF1,Max_Stress.WF2,Max_Stress.WF3,Max_Stress.WF4];

% warm med

E_WM = [E.WM1,E.WM2,E.WM3,E.WM4];
Yield_Points_WM = [Yield_Points.WM1,Yield_Points.WM2,Yield_Points.WM3,Yield_Points.WM4];
Max_Stress_WM = [Max_Stress.WM1,Max_Stress.WM2,Max_Stress.WM3,Max_Stress.WM4];

% T testing

% cold groups
T.E.E_CF_CM = ttest(E_CF,E_CM);
T.E.E_CF_CMM = ttest(E_CF,E_CMM);
T.E.E_CM_CMM = ttest(E_CM,E_CMM);
T.YP.Yield_Points_CF_CM = ttest(Yield_Points_CF,Yield_Points_CM);
T.YP.Yield_Points_CF_CMM = ttest(Yield_Points_CF,Yield_Points_CMM);
T.YP.Yield_Points_CM_CMM = ttest(Yield_Points_CM,Yield_Points_CMM);
T.MS.Max_Stress_CF_CM = ttest(Max_Stress_CF,Max_Stress_CM);
T.MS.Max_Stress_CF_CMM = ttest(Max_Stress_CF,Max_Stress_CMM);
T.MS.Max_Stress_CM_CMM = ttest(Max_Stress_CM,Max_Stress_CMM);

% room temp groups
T.E.E_RTF_RTM = ttest(E_RTF,E_RTM);
T.YP.Yield_Points_RTF_RTM = ttest(Yield_Points_RTF,Yield_Points_RTM);
T.MS.Max_Stress_RTF_RTM = ttest(Max_Stress_RTF,Max_Stress_RTM);

% warm groups
T.E.E_WF_WM = ttest2(E_WF,E_WM);
T.YP.Yield_Points_WF_WM = ttest2(Yield_Points_WF,Yield_Points_WM);
T.MS.Max_Stress_WF_WM = ttest2(Max_Stress_WF,Max_Stress_WM);

% fast groups
T.E.E_CF_RTF = ttest2(E_CF,E_RTF);
T.E.E_CF_WF = ttest2(E_CF,E_WF);
T.E.E_RTF_WF = ttest2(E_RTF,E_WF);
T.YP.Yield_Points_CF_RTF = ttest2(Yield_Points_CF,Yield_Points_RTF);
T.YP.Yield_Points_CF_WF = ttest2(Yield_Points_CF,Yield_Points_WF);
T.YP.Yield_Points_RTF_WF = ttest2(Yield_Points_RTF,Yield_Points_WF);
T.MS.Max_Stress_CF_RTF = ttest2(Max_Stress_CF,Max_Stress_RTF);
T.MS.Max_Stress_CF_WF = ttest2(Max_Stress_CF,Max_Stress_WF);
T.MS.Max_Stress_RTF_WF = ttest2(Max_Stress_RTF,Max_Stress_WF);

% Med groups (ignoring CMM for the moment, unless we have to analyze it)
T.E.E_CM_RTM = ttest2(E_CM,E_RTM);
T.E.E_CM_WM = ttest2(E_CM,E_WM);
T.E.E_RTM_WM = ttest2(E_RTM,E_WM);
T.YP.Yield_Points_CM_RTM = ttest2(Yield_Points_CM,Yield_Points_RTM);
T.YP.Yield_Points_CM_WM = ttest2(Yield_Points_CM,Yield_Points_WM);
T.YP.Yield_Points_RTM_WM = ttest2(Yield_Points_RTM,Yield_Points_WM);
T.MS.Max_Stress_CM_RTM = ttest2(Max_Stress_CM,Max_Stress_RTM);
T.MS.Max_Stress_CM_WM = ttest2(Max_Stress_CM,Max_Stress_WM);
T.MS.Max_Stress_RTM_WM = ttest2(Max_Stress_RTM,Max_Stress_WM);


%% Export Data

writetable(struct2table(T.YP), 'T_Test_YP_Data.xlsx')
writetable(struct2table(T.MS), 'T_Test_MS_Data.xlsx')
writetable(struct2table(T.E), 'T_Test_E_Data.xlsx')
writetable(struct2table(E), 'E_Data.xlsx')





function dot()
fprintf('.'); %update
end