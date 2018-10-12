function varargout = GREA(varargin)
% GREA MATLAB code for GREA.fig
%      GREA, by itself, creates a new GREA or raises the existing
%      singleton*.
%
%      H = GREA returns the handle to a new GREA or the handle to
%      the existing singleton*.
%
%      GREA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GREA.M with the given input arguments.
%
%      GREA('Property','Value',...) creates a new GREA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GREA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GREA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GREA

% Last Modified by GUIDE v2.5 06-Dec-2017 22:27:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GREA_OpeningFcn, ...
                   'gui_OutputFcn',  @GREA_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GREA is made visible.
function GREA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GREA (see VARARGIN)

clc
direct0 = [];%'D:\Dropbox\Project\mfile\Proj1_ms_software\';
set(handles.edit1,'String',[direct0,'ibrutinib_nintedanib.csv']);
axes(handles.axes1);
f = imread('step1.jpg');
imshow(f);
axes(handles.axes2);
f = imread('step2.jpg');
imshow(f);
axes(handles.axes3);
f = imread('step3.jpg');
imshow(f);
axes(handles.axes4);
f = imread('step4.jpg');
imshow(f);

% Choose default command line output for GREA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GREA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GREA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

file = get(handles.edit1,'String'); %m/z directory
file0 = strrep(file,'.csv','');
data = csvread(file); %read data file
data = sortrows(data,[2 1]); %sort data file
spl = strsplit(file0,'_'); drug1 = spl{1}; drug2 = spl{2}; %split the file name to get the names of the two drugs

trim = 0.1;
ft = 9;
lw = 3;
custom_label = 1;
mksz = 0.03; %marker size 

ndose1 = length(unique(data(:,1))); %number of doses for drug 1
ndose2 = length(unique(data(:,2))); %number of doses for drug 2

%single drug response curve for drug 1
[dose1,id1] = sort(unique(unique(data(:,1)))); %concentrations of drug 1 
dose1 = dose1(2:end); %truncated concentrations of drug 1 to exclude 0
base1 = dose1(1)^2/dose1(3); 
%on logarithmic scale, there can't be zeros. Here we set base1 to be zero,
%which is way below the lowest experiment concentration on this axis. 
surv1 = data((data(:,2) == 0),3)/100; %unaffected fraction of drug 1
surv1 = surv1(id1);  %sort surv1
surv1 = surv1(2:end); %truncated unaffected fraction of drug 1 to exclude 0

%single drug response curve for drug 2
[dose2,id2] = sort(unique(unique(data(:,2)))); %concentrations of drug 2 
dose2 = dose2(2:end); %truncated concentrations of drug 2 to exclude 0
base2 = dose2(1)^2/dose2(3);
%on logarithmic scale, there can't be zeros. Here we set base1 to be zero,
%which is way below the lowest experiment concentration on this axis. 
surv2 = data((data(:,1) == 0),3)/100; %unaffected fraction of drug 2
surv2 = surv2(id2);  %sort surv2
surv2 = surv2(2:end); %truncated unaffected fraction of drug 2 to exclude 0

%combination data for drug 1 and 2
surv12 = zeros(length(dose1),length(dose2)); %unaffected fraction in a matrix
surv12_ar = zeros(length(dose1)*length(dose2),3); %unaffected fraction in an array
k = 0; %index
for i = 1:length(dose1)
    for j = 1:length(dose2)
        k = k + 1;
        surv12(i,j) = data(intersect(find(data(:,1) == dose1(i)),...
            find(data(:,2) == dose2(j))),3)/100; %unaffected fraction in a matrix
        surv12_ar(k,:) = [dose1(i) dose2(j) surv12(i,j)]; %unaffected fraction in an array
    end 
end

%% 2 calculate Hill parameters for individual drugs
%first regression 
bot = min([surv1' surv2']); %initial guess of S0, assay background
c1 = package_envelope_hill([median(dose1) 2],[1;surv1],[base1;dose1],bot); %regression for drug 1
lam1 = c1(1); h1 = c1(2); %EC50 and Hill slope for drug 1
c2 = package_envelope_hill([median(dose2) 2],[1;surv2],[base2;dose2],bot); %regression for drug 2
lam2 = c2(1); h2 = c2(2); %EC50 and Hill slope for drug 2

%third regression with adjusted S0, assay background
c = package_envelope_s0([lam1,h1,bot,lam2,h2],surv1,dose1,surv2,dose2); %regression for S0
bot = c(3); %new S0, assay background
lam1 = c(1); h1 = c(2); %EC50 and Hill slope for drug 1
lam2 = c(4); h2 = c(5); %EC50 and Hill slope for drug 2

%% 3 visualization
%initialize figure object
%figure('position',[50 50 1500 310]);


%---------------------------  PANEL 1  ------------------------------------
%This panel provides visualization of the data and the envelope in 3-D
cla(handles.axes1)
axes(handles.axes1);

%envelope plotter
package_envelope_plotter(lam1,lam2,1,bot,h1,h2,log10(base1)-0.3,...
    log10(max(dose1))+0.3,log10(base2)-0.3,log10(max(dose2))+0.3,200,ft); 

%define 3-d markers
asp1 = (log10(max(dose1))-log10(base1)+2*trim)/1.2; %define aspect ratio for axis 1
asp2 = (log10(max(dose2))-log10(base2)+2*trim)/1.2; %define aspect ratio for axis 2
alpha = 0.7; %transparency

%single drug response for drug 1
for i = 1:length(dose1)
        hold on
        package_envelope_plotspheres(log10(dose1(i)),log10(base2),surv1(i),mksz,alpha,asp1,asp2,[0 0 0]);
end

%single drug response for drug 2
for j = 1:length(dose2)
        hold on
        package_envelope_plotspheres(log10(base1),log10(dose2(j)),surv2(j),mksz,alpha,asp1,asp2,[0 0 0]);
end

%combination data for drug 1 and 2
for k = 1:(length(dose1)*length(dose2))
    [f1,f2] = package_loewe(surv12_ar(k,1),surv12_ar(k,2),1,bot,lam1,lam2,h1,h2); %the effect calculated from generalized Loewe
    f3 = package_bliss(surv12_ar(k,1),surv12_ar(k,2),1,bot,lam1,lam2,h1,h2); %the effect calculated from Bliss
    if surv12_ar(k,3) < min([f1 f2 f3]) %red for synergy
        hold on
        package_envelope_plotspheres(log10(surv12_ar(k,1)),log10(surv12_ar(k,2)),surv12_ar(k,3),mksz,alpha,asp1,asp2,[1 0 0]);
    elseif surv12_ar(k,3) > max([f1 f2 f3]) %blue for antagonism
        hold on
        package_envelope_plotspheres(log10(surv12_ar(k,1)),log10(surv12_ar(k,2)),surv12_ar(k,3),mksz,alpha,asp1,asp2,[0 0 1]);
    else %white for additivity
        hold on 
        package_envelope_plotspheres(log10(surv12_ar(k,1)),log10(surv12_ar(k,2)),surv12_ar(k,3),mksz,alpha,asp1,asp2,[1 1 1]);
    end
end

%notations
xlim([log10(base1)-mksz*asp1 log10(max(dose1))+mksz*asp1]); %x axis limit
ylim([log10(base2)-mksz*asp2 log10(max(dose2))+mksz*asp2]); %y axis limit
xlabel([drug1,' (\muM)'],'fontsize',ft); %x axis label
ylabel([drug2,' (\muM)'],'fontsize',ft); %y axis label
zlim([0 1.2]); %z axis limit
camlight right; %define camera light
lighting phong; %define light type
%title([drug1,' and ',drug2],'fontsize',ft); %title of the figure

%replace logarithmic tick labels with algebraic ones
xlb = 2.^(round(log2(base1)):2:log2(max(dose1)));
if length(xlb)<= 3
    xlb = 2.^(round(log2(base1)):1:log2(max(dose1)));
end
ylb = 2.^(round(log2(base2)):2:log2(max(dose2)));
if length(ylb)<= 3
    ylb = 2.^(round(log2(base2)):1:log2(max(dose2)));
end
if custom_label == 1
    set(gca,'xtick',(log10(xlb)),'xticklabel',num2cell((xlb)),'ytick',(log10(ylb)),'yticklabel',num2cell((ylb)));
else
    xlabel([drug1,' (lg(\muM))'],'fontsize',ft); %x axis label
    ylabel([drug2,' (lg(\muM))'],'fontsize',ft); %y axis label
end
set(gca,'fontsize',ft);
view(150,20); %camera angle

%---------------------------  PANEL 2  ------------------------------------
%This panel provides 2-D projection of the combination data shown in panel
%1, thus no comments are shown for similar scripts. 
cla(handles.axes2)
axes(handles.axes2);

mksz = 12; %redefine marker size for this panel
%label the data points
ff = zeros(1,length(dose1)*length(dose2)); %synergy labels
gg = zeros(1,length(dose1)*length(dose2)); %antagonism labels
for k = 1:(length(dose1)*length(dose2))
    [f1,f2] = package_loewe(surv12_ar(k,1),surv12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    f3 = package_bliss(surv12_ar(k,1),surv12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    if surv12_ar(k,3) < min([f1 f2 f3])
        ff(k) = 1;
    end
    if surv12_ar(k,3) > max([f1 f2 f3])
        gg(k) = 1;
    end
end

%visualize the island
plot(log10(surv12_ar((ff==1),1)),log10(surv12_ar((ff==1),2)),'rs','markersize',mksz,'markerfacecolor','r');
hold on
plot(log10(surv12_ar((gg==1),1)),log10(surv12_ar((gg==1),2)),'bs','markersize',mksz,'markerfacecolor','b');
hold off
xlim([log10(min(dose1))-trim log10(max(dose1))+trim]);
ylim([log10(min(dose2))-trim log10(max(dose2))+trim]);
xlabel([drug1,' (\muM)'],'fontsize',ft); %x axis label
ylabel([drug2,' (\muM)'],'fontsize',ft); %y axis label

if custom_label == 1
    set(gca,'xtick',(log10(xlb)),'xticklabel',num2cell((xlb)),'ytick',(log10(ylb)),'yticklabel',num2cell((ylb)));
else
    xlabel([drug1,' (lg(\muM))'],'fontsize',ft); %x axis label
    ylabel([drug2,' (lg(\muM))'],'fontsize',ft); %y axis label
end
set(gca,'fontsize',ft);
box on

%---------------------------  PANEL 3  ------------------------------------
%This panel identifies the largest islands of synergy and antagonism using 
%connected-component labeling.
cla(handles.axes3)
axes(handles.axes3);

%label the data points
ff = zeros(1,length(dose1)*length(dose2)); %synergy labels
gg = zeros(1,length(dose1)*length(dose2)); %antagonism labels
for k = 1:(length(dose1)*length(dose2))
    [f1,f2] = package_loewe(surv12_ar(k,1),surv12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    f3 = package_bliss(surv12_ar(k,1),surv12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    if surv12_ar(k,3) < min([f1 f2 f3])
        ff(k) = 1;
    end
    if surv12_ar(k,3) > max([f1 f2 f3])
        gg(k) = 1;
    end
end

%find the largest island of synergy
ff = reshape(ff,length(dose1),length(dose2)); %reshape synergy label array to matrix
CC = bwconncomp(ff,4); %connected-component labeling
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,id] = max(numPixels);
if isempty(id) == 0
    reg = CC.PixelIdxList{id}'; %data list corresponding to the largest island of antagonism
else
    reg = [];
end

%find the largest island of antagonism
gg = reshape(gg,length(dose1),length(dose2)); %reshape antagonism label array to matrix
CC = bwconncomp(gg,4); %connected-component labeling
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,id] = max(numPixels);
if isempty(id) == 0
    reg2 = CC.PixelIdxList{id}'; %data list corresponding to the largest island of antagonism
else
    reg2 = [];
end

%examine if additivity region is larger than the other two
adv = 1 - gg - ff; %additivity region
adv = reshape(adv,length(dose1),length(dose2)); %reshape antagonism label array to matrix
CC = bwconncomp(adv,4); %connected-component labeling
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,id] = max(numPixels);
if isempty(id) == 0
    reg3 = CC.PixelIdxList{id}'; %data list corresponding to the largest island of antagonism
else
    reg3 = [];
end

%visualize the island
plot(log10(surv12_ar(reg,1)),log10(surv12_ar(reg,2)),'rs','markersize',mksz,'markerfacecolor','r');
hold on
plot(log10(surv12_ar(reg2,1)),log10(surv12_ar(reg2,2)),'bs','markersize',mksz,'markerfacecolor','b');
hold off
xlim([log10(min(dose1))-trim log10(max(dose1))+trim]);
ylim([log10(min(dose2))-trim log10(max(dose2))+trim]);
xlabel([drug1,' (\muM)'],'fontsize',ft); %x axis label
ylabel([drug2,' (\muM)'],'fontsize',ft); %y axis label

if custom_label == 1
    set(gca,'xtick',(log10(xlb)),'xticklabel',num2cell((xlb)),'ytick',(log10(ylb)),'yticklabel',num2cell((ylb)));
else
    xlabel([drug1,' (lg(\muM))'],'fontsize',ft); %x axis label
    ylabel([drug2,' (lg(\muM))'],'fontsize',ft); %y axis label
end
set(gca,'fontsize',ft);

%---------------------------  PANEL 4  ------------------------------------
%This panel visualizes the difference of the measured data from the 
%envelope for synergy and antagonism islands and calculates the synergy 
%index (SI) and antagonism index (AI). 
cla(handles.axes4)
axes(handles.axes4);

%added 6/29/2018, linear scale correction
d1p = [dose1(1)^2/dose1(2);dose1];
d2p = [dose2(1)^2/dose2(2);dose2];

%envelope plotter
mksz = 0.015; %redefine marker size for this panel
package_envelope_plotter(lam1,lam2,1,bot,h1,h2,log10(base1)-0.3,...
    log10(max(dose1))+0.3,log10(base2)-0.3,log10(max(dose2))+0.3,200,ft);

hold on
%difference of the measured data from the envelope for synergy
%reg = setdiff(reg,idss);
%reg2 = setdiff(reg2,idss);

nn = 0;
area1 = 0;
if isempty(reg) == 0
dif = zeros(1,length(reg));
for k = reg
    nn = nn + 1;
    package_envelope_plotspheres(log10(surv12_ar(k,1)),log10(surv12_ar(k,2)),surv12_ar(k,3),mksz,alpha,asp1,asp2,[1 0 0]);
    hold on
    [f1,f2] = package_loewe(surv12_ar(k,1),surv12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    f3 = package_bliss(surv12_ar(k,1),surv12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    flow = min([f1 f2 f3]);
    plot3([log10(surv12_ar(k,1)) log10(surv12_ar(k,1))],[log10(surv12_ar(k,2)) log10(surv12_ar(k,2))],...
        [surv12_ar(k,3) flow],'linewidth',lw,'color','r');
    hold on
    package_envelope_plotspheres(log10(surv12_ar(k,1)),log10(surv12_ar(k,2)),flow,mksz,alpha,asp1,asp2,[1 0 0]);
    hold on
    id1 = find(d1p == surv12_ar(k,1));
    id2 = find(d2p == surv12_ar(k,2));
    dif(nn) = abs(surv12_ar(k,3) - flow)*log(d1p(id1)/d1p(id1-1))*log(d2p(id2)/d2p(id2-1));
    area1 = area1 + log(d1p(id1)/d1p(id1-1))*log(d2p(id2)/d2p(id2-1));
end
else 
dif = 0;   
end

%difference of the measured data from the envelope for antagonism
nn = 0;
area2 = 0;
if isempty(reg2) == 0
dif2= zeros(1,length(reg2));
for k = reg2
    nn = nn + 1;
    package_envelope_plotspheres(log10(surv12_ar(k,1)),log10(surv12_ar(k,2)),surv12_ar(k,3),mksz,alpha,asp1,asp2,[0 0 1]);
    hold on
    [f1,f2] = package_loewe(surv12_ar(k,1),surv12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    f3 = package_bliss(surv12_ar(k,1),surv12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    flow = max([f1 f2 f3]);
    plot3([log10(surv12_ar(k,1)) log10(surv12_ar(k,1))],[log10(surv12_ar(k,2)) log10(surv12_ar(k,2))],...
        [surv12_ar(k,3) flow],'linewidth',lw,'color','b');
    hold on
    package_envelope_plotspheres(log10(surv12_ar(k,1)),log10(surv12_ar(k,2)),flow,mksz,alpha,asp1,asp2,[0 0 1]);
    hold on
    id1 = find(d1p == surv12_ar(k,1));
    id2 = find(d2p == surv12_ar(k,2));
    dif2(nn) = abs(surv12_ar(k,3) - flow) * log(d1p(id1)/d1p(id1-1))*log(d2p(id2)/d2p(id2-1));
    area2 = area2 + log(d1p(id1)/d1p(id1-1))*log(d2p(id2)/d2p(id2-1));
end
else 
dif2 = 0;  
end

%difference of the measured data from the envelope for antagonism
nn = 0;
area3 = 0;
if isempty(reg3) == 0
for k = reg3
    nn = nn + 1;
    id1 = find(d1p == surv12_ar(k,1));
    id2 = find(d2p == surv12_ar(k,2));
    area3 = area3 + log(d1p(id1)/d1p(id1-1))*log(d2p(id2)/d2p(id2-1));
end
end
hold off

if area3 >= max(area1,area2)
    dif = 0;
    dif2 = 0;
end

tarea = 0;
for i = 2:length(d1p)
    for j = 2:length(d2p)
        tarea = tarea + log(d1p(i)/d1p(i-1))*log(d2p(j)/d2p(j-1));
    end
end
si = sum(dif)/tarea; %SI
ai = sum(dif2)/tarea; %AI

set(handles.edit2,'String',num2str(si,3));
set(handles.edit3,'String',num2str(ai,3));

xlim([log10(base1)-mksz*asp1 log10(max(dose1))+mksz*asp1]);
ylim([log10(base2)-mksz*asp2 log10(max(dose2))+mksz*asp2]);
xlabel([drug1,' (\muM)'],'fontsize',ft); %x axis label
ylabel([drug2,' (\muM)'],'fontsize',ft); %y axis label
zlim([0 1.2]);
camlight right; 
lighting phong;
%title(['SI = ',num2str(si,2),', AI = ',num2str(ai,2)],'fontsize',ft);

if custom_label == 1
    set(gca,'xtick',(log10(xlb)),'xticklabel',num2cell((xlb)),'ytick',(log10(ylb)),'yticklabel',num2cell((ylb)));
else
    xlabel([drug1,' (lg(\muM))'],'fontsize',ft); %x axis label
    ylabel([drug2,' (lg(\muM))'],'fontsize',ft); %y axis label
end
set(gca,'fontsize',ft);
view(150,20);

[alpha,beta,gamma] = RSRA_package(data,0,0.5); %RSRA
set(handles.edit4,'String',num2str(alpha,3));
set(handles.edit5,'String',num2str(beta,3));
set(handles.edit6,'String',num2str(gamma,3));
