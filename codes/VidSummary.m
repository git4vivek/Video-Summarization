function varargout = VidSummary(varargin)
% VIDSUMMARY M-file for VidSummary.fig
%      VIDSUMMARY, by itself, creates a new VIDSUMMARY or raises the existing
%      singleton*.
%
%      H = VIDSUMMARY returns the handle to a new VIDSUMMARY or the handle to
%      the existing singleton*.
%
%      VIDSUMMARY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIDSUMMARY.M with the given input arguments.
%
%      VIDSUMMARY('Property','Value',...) creates a new VIDSUMMARY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before VidSummary_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to VidSummary_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help VidSummary

% Last Modified by GUIDE v2.5 15-May-2012 09:58:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @VidSummary_OpeningFcn, ...
                   'gui_OutputFcn',  @VidSummary_OutputFcn, ...
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


% --- Executes just before VidSummary is made visible.
function VidSummary_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VidSummary (see VARARGIN)

% Choose default command line output for VidSummary
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
handles.FileName = '';
handles.extension = '';
handles.KeyFrame = 0;
handles.Mov = 0;
handles.FolderName = '';
handles.Obj= 0;
handles.Omov=0;
% UIWAIT makes VidSummary wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VidSummary_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Browse.
function Browse_Callback(hObject, eventdata, handles)
% hObject    handle to Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.FileName = uigetfile('*.avi');
    guidata(hObject,handles);
    
    


% --- Executes on button press in Summarize.
function Summarize_Callback(hObject, eventdata, handles)
% hObject    handle to Summarize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


obj = mmreader(handles.FileName);
handles.Omov=obj;
nFrames = obj.NumberOfFrames;
p=1;

for i= 1 : fix(nFrames/300)+1
   if i==1
       o=5;
       r=300;
   else
       o= (5+(i-1)*300)-100;
       r= min(i*300, nFrames)-300;
   end
   if i == (fix(nFrames/300)+1)
       r = min(i*300, nFrames);
   end
   
for k = o : 5 : r
    img = read(obj,k);
    %img = rgb2gray(img);
    %img = rgb2hsv(img);
    img = imresize(img,0.05);
    img = im2double(img);
    %img = histeq(img);
    [m n] = size(img);
    img = reshape(img,1,m*n);
    
    for j = 1:m*n
            c1(((k-o)/5)+1,j) = img(1,j);
    end
   
end

key = myfun(c1);
[r c] = size(key);
for l=1:c
keyframe(p) = o-5+key(l)*5;
p=p+1;
end
c1 = [];
i
end
keyframe=unique(keyframe)

[r c] = size(keyframe);
nFrames = c;
vidHeight = obj.Height;
vidWidth = obj.Width;

% Preallocate movie structure.
mov(1:nFrames) = struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),'colormap', []);

% Read one frame at a time.
for k = 1 : nFrames
    mov(k).cdata = read(obj, keyframe(k));
end
handles.Mov = mov;
cd('C:\Users\Vivek\Desktop\Video-Summarization\summaries');
FolderName = ['SumVid' handles.FileName]; 
movie2avi(mov, FolderName);
handles.FolderName = FolderName;
cd('C:\Users\Vivek\Desktop\Video-Summarization\codes');
beep
beep
beep
handles.KeyFrame = keyframe;
handles.Obj = obj;
guidata(hObject,handles);






% --- Executes on button press in Original.
function Original_Callback(hObject, eventdata, handles)
% hObject    handle to Original (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

implay(handles.FileName);
% obj = handles.Omov;
% nFrames = obj.NumberOfFrames;
% vidHeight = handles.Omov.Height;
% vidWidth = handles.Omov.Width;
% mov(1:nFrames) = struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),'colormap', []);
% 
% % Read one frame at a time.
% for k = 1 : nFrames
%     mov(k).cdata = read(obj, k);
% end
% % Size a figure based on the video's width and height.
% 
% % Play back the movie once at the video's frame rate.
% hf = figure;
% set(hf, 'position', [150 150 vidWidth vidHeight])
% 
% movie(hf,mov);


% --- Executes on button press in VidSummary.
function VidSummary_Callback(hObject, eventdata, handles)
% hObject    handle to VidSummary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd('C:\Users\Vivek\Desktop\Video-Summarization\summaries');
sumname = ['SumVid' handles.FileName];
implay(sumname);
cd('C:\Users\Vivek\Desktop\Video-Summarization\codes');
% mov = handles.Mov;
% vidHeight = handles.Obj.Height;
% vidWidth = handles.Obj.Width;
% 
% % Size a figure based on the video's width and height.
% 
% % Play back the movie once at the video's frame rate.
% hf = figure;
% set(hf, 'position', [150 150 vidWidth vidHeight])
% 
% movie(hf,mov, 1, 8);





% --- Executes on button press in ImageSummary.
function ImageSummary_Callback(hObject, eventdata, handles)
% hObject    handle to ImageSummary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%obj = handles.Mov;
modifiedStr = strrep(handles.FileName,'.avi','');
obj = mmreader(handles.FileName);
    %vid = read(obj);
    %nFrame = size(vid,4);  
ImageSum = ['SumImage' modifiedStr];    
    [n m] = size(handles.KeyFrame);
    handles.KeyFrame
    cd ('C:\Users\Vivek\Desktop\Video-Summarization\summaries');
    mkdir(ImageSum);
    cd(ImageSum);
    for k= 1 : 9
        newname = strcat(num2str(handles.KeyFrame(k)),'.jpg');
        imwrite(obj.read(handles.KeyFrame(min((fix((m+1)/18)+(k-1)*(fix(m/9)+1)),m))),newname);
        %imwrite(vid(:,:,:,(handles.KeyFrame(k))), newname);  
    end
   
fileFolder = fullfile('C:','Users','Vivek','Desktop','Video-Summarization','summaries',ImageSum);
dirOutput = dir(fullfile(fileFolder,'*.jpg'));
fileNames = {dirOutput.name}'
montage(fileNames);
 cd ('C:\Users\Vivek\Desktop\Video-Summarization\codes');




% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nFrames = handles.Obj.NumberOfFrames;
[n m]= size(handles.KeyFrame);

Duration = handles.Obj.Duration;
Framerate = handles.Obj.FrameRate;
MovieDur = m/8;
f = figure('Position',[150 150 500 100]);
rnames = {'Original Video','Summarized Video'};
cnames = {'No. Of Frames','Frame Rate(per sec)','Duration(in sec)'};
dat =  {nFrames , Framerate, Duration; m , 5 , MovieDur; };
t = uitable('Data',dat,'ColumnName',cnames, 'RowName',rnames,'Position',[20 20 460 80]);
