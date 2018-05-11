function varargout = gui(varargin)
clc
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 16-Oct-2016 18:23:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
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


% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;

handles.flagMovie=0;
% Update handles structure
guidata(hObject, handles);

set(handles.pbuttonRun,'BackgroundColor','g')


% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function etextCFL_Callback(hObject, eventdata, handles)
% hObject    handle to etextCFL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etextCFL as text
%        str2double(get(hObject,'String')) returns contents of etextCFL as a double
flag = handles.IC;
if handles.IC == 1 
    sodConditions(handles);
else
    shuConditions(handles);
end


% --- Executes during object creation, after setting all properties.
function etextCFL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etextCFL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etextDeltax_Callback(hObject, eventdata, handles)
% hObject    handle to etextDeltax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etextDeltax as text
%        str2double(get(hObject,'String')) returns contents of etextDeltax as a double
flag = handles.IC;
if handles.IC == 1 
    sodConditions(handles);
else
    shuConditions(handles);
end


% --- Executes during object creation, after setting all properties.
function etextDeltax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etextDeltax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etextDeltat_Callback(hObject, eventdata, handles)
% hObject    handle to etextDeltat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etextDeltat as text
%        str2double(get(hObject,'String')) returns contents of etextDeltat as a double



% --- Executes during object creation, after setting all properties.
function etextDeltat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etextDeltat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkRoe.
function checkRoe_Callback(hObject, eventdata, handles)
% hObject    handle to checkRoe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkRoe


% --- Executes on button press in checkWeno.
function checkWeno_Callback(hObject, eventdata, handles)
% hObject    handle to checkWeno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkWeno


% --- Executes on button press in checkRoeWeno.
function checkRoeWeno_Callback(hObject, eventdata, handles)
% hObject    handle to checkRoeWeno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkRoeWeno


% --- Executes on button press in checkCompact.
function checkCompact_Callback(hObject, eventdata, handles)
% hObject    handle to checkCompact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkCompact


% --- Executes on button press in pbuttonRun.
function pbuttonRun_Callback(hObject, eventdata, handles)
% hObject    handle to pbuttonRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pbuttonPlay,'BackgroundColor',[0.941, 0.941, 0.941]);
set(handles.pbuttonExpAni,'BackgroundColor',[0.941, 0.941, 0.941]);
set(handles.pbuttonRun,'BackgroundColor','y')
ratio = str2num(get(handles.etextRatio,'String'));

inp.IC = handles.IC;
CFL = str2num(get(handles.etextCFL,'String'));
inp.DELTA_X = str2num(get(handles.etextDeltax,'String'));

flagDone=0; % Flag for calculating the analytic solution of Sod 
if (get(handles.checkRoe,'Value')) == 1 
   fComp=0;
   fRoe=1;
   inp.WENO=0;
   codeGui
end

if (get(handles.checkWeno,'Value')) == 1 
   fComp=0;
   fRoe=0;
   inp.WENO=1;
   codeGui
end

if (get(handles.checkRoeWeno,'Value')) == 1 
   fComp=0;
   fRoe=1;
   inp.WENO=1;
   codeGui
end

if (get(handles.checkCompact,'Value')) == 1 
   fComp=1;
   fRoe=0;
   inp.WENO=0;
   if get(handles.checkVisc,'Value')==1
       inp.SHEAR = str2num(get(handles.etextVisc,'String'));
   else
       inp.SHEAR=0;
   end
   codeGui
end


if flagOk == 1
    set(handles.pbuttonRun,'BackgroundColor','g')
else
    set(handles.pbuttonRun,'BackgroundColor','r')
    beep
end


% --- Executes on button press in pbuttonImport.
function pbuttonImport_Callback(hObject, eventdata, handles)
% hObject    handle to pbuttonImport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in .
function pbuttonPlay_Callback(hObject, eventdata, handles)
% hObject    handle to pbuttonPlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.pbuttonRun,'BackgroundColor',[0.941, 0.941, 0.941]);
if (handles.flagMovie==0)
   set(handles.pbuttonExpAni,'BackgroundColor',[0.941, 0.941, 0.941]); 
end
set(handles.pbuttonPlay,'BackgroundColor','y');
set(handles.stextImpStat,'Visible','off');

rat = str2num(get(handles.etextRatio,'String'));
cfl = str2num(get(handles.etextCFL,'String'));
dx = str2num(get(handles.etextDeltax,'String'));

% Identifies what's the name of the files
n=0;
filename{1}=[];
if handles.IC ~= 1
    rat=0;
end
if (get(handles.checkRoe,'Value')) == 1 
   flag(1)=0; % Comp
   flag(2)=1; % Roe 
   flag(3)=0; % Weno
   n=n+1;
   filename{n} = makeFilename (flag, handles.IC, cfl, dx, 0, rat);
end

if (get(handles.checkWeno,'Value')) == 1 
   flag(1)=0;
   flag(2)=0;
   flag(3)=1;
   n=n+1;
   filename{n} = makeFilename (flag, handles.IC, cfl, dx, 0, rat);
end

if (get(handles.checkRoeWeno,'Value')) == 1 
   flag(1)=0;
   flag(2)=1;
   flag(3)=1;
   n=n+1;
   filename{n} = makeFilename (flag, handles.IC, cfl, dx, 0, rat);
end

if (get(handles.checkCompact,'Value')) == 1 
   flag(1)=1;
   flag(2)=0;
   flag(3)=0;
   if get(handles.checkVisc,'Value')==1
       vis = str2num(get(handles.etextVisc,'String'));
   else
       vis=0;
   end
   n=n+1;
   filename{n} = makeFilename (flag, handles.IC, cfl, dx, vis, rat);
end

% Check if files exist
flag=1;
for i =1:n
    if exist(filename{i}) == 0
        set(handles.pbuttonPlay,'BackgroundColor','r')
        set(handles.stextImpStat,'String','Datafile not found')
        set(handles.stextImpStat,'Visible','on')
        beep
        flag=0;
    end
end

% Load data
if flag==1
    set(handles.stextImpStat,'String','Datafiles loaded')    
    set(handles.stextImpStat,'Visible','on')
    prop = get(handles.pmenuAxis,'Value');
    for i=1:n
        temp = load(filename{i},'U');        
        if prop == 4
            f(i,:,:) = 0.4 * (temp.U(3,:,:) - 0.5*temp.U(1,:,:).*(temp.U(2,:,:) ./ temp.U(1,:,:)).^2);
        elseif prop==5
            f(i,:,:) = temp.U(2,:,:) ./ temp.U(1,:,:);
        elseif prop==6
            f(i,:,:) = abs(temp.U(2,:,:) ./ temp.U(1,:,:)) ./ sqrt(1.4*(0.4 * (temp.U(3,:,:) - 0.5*temp.U(1,:,:).*(temp.U(2,:,:) ./ temp.U(1,:,:)).^2))./temp.U(1,:,:));
        else
            f(i,:,:) = temp.U(prop,:,:);
        end
    end
    temp = load(filename{i},'x');
    x = temp.x;
    clear temp;
    
    % Load reference solution, if it exists
    refname = strcat(filename{1}(4:6), 'Rat', get(handles.etextRatio,'String'),'Reference.mat');
    flagAnalytic=0;
    if exist(refname) ~= 0
        temp = load(refname);
        flagAnalytic = 1;
        xa = temp.x;
        if prop == 4
            fa = 0.4 * (temp.f(3,:,:) - 0.5*temp.f(1,:,:).*(temp.f(2,:,:) ./ temp.f(1,:,:)).^2);
        elseif prop==5
            fa = temp.f(2,:,:) ./ temp.f(1,:,:);
        elseif prop==6
            fa = abs(temp.f(2,:,:) ./ temp.f(1,:,:)) ./ sqrt(1.4*(0.4 * (temp.f(3,:,:) - 0.5*temp.f(1,:,:).*(temp.f(2,:,:) ./ temp.f(1,:,:)).^2))./temp.f(1,:,:));
        else
            fa = temp.f(prop,:,:);
        end
        nfa = size(fa);
        

    end

% Animation
color = 'gcym';
nf = size(f);

velcoef = round((1-get(handles.sliderVel,'Value'))/0.5*45);
if velcoef==0
    velcoef=1;
end
ptime=0.01;
if (get(handles.sliderVel,'Value')) < 0.3
    ptime = 0.02;
end
if (get(handles.sliderVel,'Value'))<0.15
    ptime=0.1;
end

if (handles.flagMovie==1)
    handles.H=figure('units','pixels','outerposition',[0 0 1280 720]);
    set(handles.H,'Color','w')
    ptime=0;
end

step=nf(3)/velcoef;

supy = max(max(f(:,:,nf(3))));
infy = min(min(f(:,:,nf(3))));

for j=1:velcoef
    hold off
    if flagAnalytic==1
       stepAna = nfa(3)/velcoef;
       plotAnalitic(handles, xa, fa(1,:,round(j*stepAna)), 'w');
       hold on 
    end
    for i=1:n
        plotVariable(handles, x, f(i,:,round(j*step)), color(i))
        if max(max(f(:,:,round(j*step))))>supy
            supy=max(max(f(:,:,round(j*step))));
        end
        
        if min(min(f(:,:,round(j*step)))) < infy
            infy = min(min(f(:,:,round(j*step))));
        end
        axis([ x(1), x(end), infy, supy ])
        hold on
    end
    if handles.flagMovie==1
        M(j)=getframe(gcf);
    end
    pause(ptime)
    
end

if flagAnalytic==0
    leg = cell(1,n);
    for i = 1:n
        leg{i} = filename{i}(1:3);
    end
else
    leg = cell(1,n+1);
    leg(1) = {'Reference'};
    for i = 2:n+1
        leg{i} = filename{i-1}(1:3);
    end
end
legend (leg)
set(legend,'Textcolor','w')
set(legend)
set(handles.graph,'FontSize',16);
M(velcoef:velcoef+round(velcoef*0.4))=getframe(gcf);

set(handles.pbuttonPlay,'BackgroundColor','g')
set(handles.stextImpStat,'Visible','off')

if handles.flagMovie==1
    % Avoid bugs of getframe
    nM = size(M(1).cdata);
    nM(1)=nM(1)-1;
    nM(2)=nM(2)-1;
    for i=1:velcoef+round(velcoef*0.4)
        temp = size(M(i).cdata);
        if (temp(1) ~= nM(1))
            fat = temp(1)-nM(1);
            M(i).cdata(nM(1)-fat+1:nM(1),:,:) = []; 
        end
        if (temp(2) ~= nM(2))
            fat = temp(2)-nM(2);
            M(i).cdata(:,nM(2)-fat+1:nM(2),:) = []; 
        end       
    end
    
    writerObj =  VideoWriter('newfile','MPEG-4');
    writerObj.Quality = 100;
    open(writerObj);
    writeVideo(writerObj,M);
    close(writerObj);
    set(handles.pbuttonExpAni,'BackgroundColor','g')
    handles.flagMovie=0;
    close (handles.H)
    guidata(handles.output,handles)
end

end

function [filename] = makeFilename (flag, ico, cfl, dx, vis, rat)

if flag(1) == 1 && flag(2) ==0 && flag(3) ==0
    name{1} = 'Com';
end

if flag(1) == 0 && flag(2) == 1 && flag(3) == 0 
    name{1} = 'Roe';
    vis=0;
end

if flag(1) == 0 && flag(2) == 0 && flag(3) == 1 
    name{1} = 'Wen';
    vis=0;
end

if flag(1) ==0 && flag(2) == 1 && flag(3) == 1 
    name{1} = 'Rwe';
    vis=0;
end

name{2} = { 'Sod', 'Shu', 'Rsh', 'Prs' };
name{2} = name{2}(ico);

name{3} =  'Cfl' ;

name{4} =  'Dx' ;

name{5} =  'Vis' ;

name{6} =  'Rat' ;

filename = strcat(name(1), ...
                  name{2}(1),...
                  name(3), num2str(cfl), ...
                  name(4), num2str(dx), ...
                  name(5), num2str(vis), ...
                  name(6), num2str(rat));
              
filename = strrep(filename, '.', '');
filename = strcat(filename{1},'.mat');


% --- Executes on button press in pbuttonExpFig.
function pbuttonExpFig_Callback(hObject, eventdata, handles)
% hObject    handle to pbuttonExpFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(0,'showhiddenhandles','on') % Make the GUI figure handle visible
h = findobj(gcf,'type','axes'); % Find the axes object in the GUI
f1 = figure('units','normalized','outerposition',[0 0 1 1]); % Open a new figure with handle f1
set(gcf,'Color','w')
s = copyobj(h,f1); % Copy axes object h into figure f1
export_fig figure.pdf



% --- Executes on selection change in pmenuAxis.
function pmenuAxis_Callback(hObject, eventdata, handles)
% hObject    handle to pmenuAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pmenuAxis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pmenuAxis

fSod = get(handles.rbuttonSod, 'Value');
fShu = get(handles.rbuttonShu, 'Value');
fRefShu = get(handles.rbuttonRefShu, 'Value');
fParRefShu = get(handles.rbuttonParRefShu,'Value');
if fSod == 1
    sodConditions(handles)
elseif fShu == 1 
    shuConditions(handles)
elseif fRefShu == 1
    refShuConditions(handles)
elseif fParRefShu == 1
    parRefShuConditions(handles)
end


% --- Executes during object creation, after setting all properties.
function pmenuAxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmenuAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbuttonSod.
function rbuttonSod_Callback(hObject, eventdata, handles)
% hObject    handle to rbuttonSod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbuttonSod
set(handles.rbuttonSod,'Value',1)
set(handles.rbuttonShu,'Value',0)
set(handles.rbuttonRefShu,'Value',0)
set(handles.rbuttonParRefShu,'Value',0)
set(handles.stextRatio,'Visible','on')
set(handles.etextRatio,'Visible','on')
set(handles.sliderVel,'Value',0.5);
set(handles.etextRatio,'String','10');


set(handles.etextDeltax,'String', 0.1)
sodConditions(handles)

handles.IC = 1;
guidata(handles.output, handles)



function sodConditions (handles)
cfl = str2num(get(handles.etextCFL,'String'));
dx = str2num(get(handles.etextDeltax,'String'));
ratio = str2num(get(handles.etextRatio,'String'));

infx = -5;
supx = 5;
Pl = 0.1*ratio;
rhol = Pl;
Co = (1.4*Pl./rhol)^0.5;
dt = cfl*dx / ( 0 + Co );
rhor=0.125;
Pr=0.1;

x = infx:dx:supx;
sizex = length(x);

U(1,1:round(sizex/2))=rhol;
U(1,round(sizex/2):sizex)=rhor;
U(2,:) = 0;
Po(1:round(sizex/2))=Pl;
Po(round(sizex/2):sizex)=Pr;
U(3,:)=Po ./ (1.4-1);
U(4,:) = Po;
U(5,:) = U(2,:)./U(1,:);
U(6,:) = abs(U(5,:))./ (U(4,:)*1.4/U(1,:))^0.5;

set(handles.etextDeltat,'String',num2str(dt))

hold off


prop = get(handles.pmenuAxis,'Value');
plotVariable(handles, x, U(prop,:), 'w')
guidata(handles.output, handles);


% --- Executes on button press in rbuttonShu.
function rbuttonShu_Callback(hObject, eventdata, handles)
% hObject    handle to rbuttonShu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbuttonShu
set(handles.rbuttonShu,'Value',1)
set(handles.rbuttonSod,'Value',0)
set(handles.rbuttonRefShu,'Value',0)
set(handles.rbuttonParRefShu,'Value',0)
set(handles.stextRatio,'Visible','off')
set(handles.etextRatio,'Visible','off')
set(handles.etextRatio,'String','0');
set(handles.sliderVel,'Value',0.5);

set(handles.etextDeltax,'String', 0.01/5)
shuConditions(handles)
handles.IC = 2;
guidata(handles.output, handles)


function shuConditions (handles)
cfl = str2num(get(handles.etextCFL,'String'));
dx = str2num(get(handles.etextDeltax,'String'));

Co = (1.4*10.3333./3.857143)^0.5;
dt= cfl*dx / ( 2.629369 + Co );

infx = 0;
supx = 1;
x = infx:dx:supx;
sizex=length(x);

lambda = 1/8;
    
i=0;
flag=0;
while flag~=1
    i=i+1;
    if x(i)<=1/8 && x(i+1)>=1/8
        div=i;
        flag=1;
    end
end

U=zeros(6,sizex);
Po=U(1,:);
veloci=Po;

U(1,1:div-1) = 3.857143;
U(1,div:end) = 1+0.2*sin(1/lambda*2*pi*x(div:end));

veloci(1:div-1) = 2.629369;
veloci(div:sizex) = 0;

U(2,:) = U(1,:).*veloci;

Po(1:div-1) = 10.3333;
Po(div:end) = 1;

U(3,:) = Po ./ (1.4-1) + 0.5*U(2,:).*veloci;

U(4,:) = Po;

U(5,:) = U(2,:)./ U(1,:);

U(6,:) = abs(U(5,:)) ./ (U(4,:)*1.4/U(1,:))^0.5;

hold off
prop = get(handles.pmenuAxis,'Value');
plotVariable(handles, x, U(prop,:), 'w')
set(handles.etextDeltat,'String',num2str(dt))

% --- Executes on button press in rbuttonRefShu.
function rbuttonRefShu_Callback(hObject, eventdata, handles)
% hObject    handle to rbuttonRefShu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbuttonRefShu
set(handles.rbuttonRefShu,'Value',1)
set(handles.rbuttonSod,'Value',1)
set(handles.rbuttonShu,'Value',0)
set(handles.rbuttonSod,'Value',0)
set(handles.rbuttonParRefShu,'Value',0)
set(handles.stextRatio,'Visible','off')
set(handles.etextRatio,'Visible','off')
set(handles.etextRatio,'String','0');
set(handles.sliderVel,'Value',0.5);

set(handles.etextDeltax,'String', 0.01/5)
refShuConditions(handles)
handles.IC = 3;
guidata(handles.output, handles)

function refShuConditions (handles)
cfl = str2num(get(handles.etextCFL,'String'));
dx = str2num(get(handles.etextDeltax,'String'));

Co = (1.4*10.3333./3.857143)^0.5;
dt= cfl*dx / ( 2.629369 + Co );

infx = 0;
supx = 1;
x = infx:dx:supx;
sizex=length(x);

lambda = 1/8;
    
i=0;
flag=0;
while flag~=1
    i=i+1;
    if x(i)<=1/8 && x(i+1)>=1/8
        div=i;
        flag=1;
    end
end

U=zeros(6,sizex);
Po=U(1,:);
veloci=Po;

U(1,1:div-1) = 3.857143;
U(1,sizex-div:sizex) = 3.857143/2*2;
U(1,div:sizex-div-1) = 1+0.2*sin(1/lambda*2*pi*x(div:sizex-div-1));

veloci(:) = 0;
veloci(1:div-1) = 2.629369;
veloci(sizex-div:sizex) = -2.629369;

U(2,:) = U(1,:).*veloci;

Po(:) = 1;
Po(1:div-1) = 10.3333;
Po(sizex-div:sizex) = 10.3333/2*2;

U(3,:,1) = Po ./ (1.4-1) + 0.5*U(2,:).*veloci;

U(4,:) = Po;

U(5,:) = U(2,:)./ U(1,:);

U(6,:) = abs(U(5,:)) ./ (U(4,:)*1.4/U(1,:))^0.5;

hold off
prop = get(handles.pmenuAxis,'Value');
plotVariable(handles, x, U(prop,:), 'w')
set(handles.etextDeltat,'String',num2str(dt))


% --- Executes on button press in rbuttonParRefShu.
function rbuttonParRefShu_Callback(hObject, eventdata, handles)
% hObject    handle to rbuttonParRefShu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbuttonParRefShu
set(handles.rbuttonParRefShu,'Value',1)
set(handles.rbuttonSod,'Value',1)
set(handles.rbuttonShu,'Value',0)
set(handles.rbuttonRefShu,'Value',0)
set(handles.rbuttonSod,'Value',0)
set(handles.stextRatio,'Visible','off')
set(handles.etextRatio,'Visible','off')
set(handles.etextRatio,'String','0');
set(handles.sliderVel,'Value',0.5);

set(handles.etextDeltax,'String', 0.01/5)
parRefShuConditions(handles)
handles.IC = 4;
guidata(handles.output, handles)


function parRefShuConditions (handles)
cfl = str2num(get(handles.etextCFL,'String'));
dx = str2num(get(handles.etextDeltax,'String'));

Co = (1.4*10.3333./3.857143)^0.5;
dt= cfl*dx / ( 2.629369 + Co );

infx = 0;
supx = 1;
x = infx:dx:supx;
sizex=length(x);

lambda = 1/8;
    
i=0;
flag=0;
while flag~=1
    i=i+1;
    if x(i)<=1/8 && x(i+1)>=1/8
        div=i;
        flag=1;
    end
end

U=zeros(6,sizex);
Po=U(1,:);
veloci=Po;

U(1,1:div-1) = 3.857143;
U(1,sizex-div:sizex) = 3.857143/2;
U(1,div:sizex-div-1) = 1+0.2*sin(1/lambda*2*pi*x(div:sizex-div-1));

veloci(:) = 0;
veloci(1:div-1) = 2.629369;
veloci(sizex-div:sizex) = -2.629369/2;

U(2,:) = U(1,:).*veloci;

Po(:) = 1;
Po(1:div-1) = 10.3333;
Po(sizex-div:sizex) = 10.3333/2;

U(3,:,1) = Po ./ (1.4-1) + 0.5*U(2,:).*veloci;

U(4,:) = Po;

U(5,:) = U(2,:)./ U(1,:);

U(6,:) = abs(U(5,:)) ./ (U(4,:)*1.4/U(1,:))^0.5;

hold off
prop = get(handles.pmenuAxis,'Value');
plotVariable(handles, x, U(prop,:), 'w')
set(handles.etextDeltat,'String',num2str(dt))


% --- Executes on button press in checkVisc.
function checkVisc_Callback(hObject, eventdata, handles)
% hObject    handle to checkVisc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkVisc



function etextVisc_Callback(hObject, eventdata, handles)
% hObject    handle to etextVisc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etextVisc as text
%        str2double(get(hObject,'String')) returns contents of etextVisc as a double


% --- Executes during object creation, after setting all properties.
function etextVisc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etextVisc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etextRatio_Callback(hObject, eventdata, handles)
% hObject    handle to etextRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etextRatio as text
%        str2double(get(hObject,'String')) returns contents of etextRatio as a double
sodConditions(handles);

% --- Executes during object creation, after setting all properties.
function etextRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etextRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotVariable (handles, x, U, color)

i = get(handles.pmenuAxis,'Value');
plot(x, U,strcat(color,'-'),'MarkerFaceColor',color,'MarkerSize',2,'LineWidth',2)
grid on
% axis([ x(1), x(end), min(U), max(U)])
set(gca,'Xcolor',[0.3 0.3 0.3]);
set(gca,'Ycolor',[0.3 0.3 0.3]);
set(gca,'GridLineStyle','-')
set(gca,'Color','k')
xlabel('x','FontSize',16)

ystring = get(handles.pmenuAxis,'String');
ylabel(ystring{i},'FontSize',16)

function plotAnalitic (handles, x, U, color)

i = get(handles.pmenuAxis,'Value');
plot(x, U,strcat(color,'-'),'LineWidth',1.5)
grid on
% axis([ x(1), x(end), min(U), max(U)])
set(gca,'Xcolor',[0.3 0.3 0.3]);
set(gca,'Ycolor',[0.3 0.3 0.3]);
set(gca,'GridLineStyle','-')
set(gca,'Color','k')



% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function sliderVel_Callback(hObject, eventdata, handles)
% hObject    handle to sliderVel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderVel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderVel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pbuttonExpAni.
function pbuttonExpAni_Callback(hObject, eventdata, handles)
% hObject    handle to pbuttonExpAni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pbuttonRun,'BackgroundColor',[0.941, 0.941, 0.941]);
set(handles.pbuttonExpAni,'BackgroundColor','y')
handles.flagMovie=1;
guidata(handles.output,handles)
pbuttonPlay_Callback(handles.output, 0, handles)
