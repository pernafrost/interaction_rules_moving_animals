function varargout = matlab_visualization(varargin)
% function matlab_visualization
% The function opens a visualization window to play with the interaction
% rules of flocking particles.
% The focal individual aims at reaching a position on the side of a close
% neighbour (or in front of the neighobur, or other, depending on the
% parameters).
%
% Written by:
% Andrea Perna
% http://www.perna.fr
%
% Date:
% 2014 / 04 / 18



% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @matlab_visualization_OpeningFcn, ...
    'gui_OutputFcn',  @matlab_visualization_OutputFcn, ...
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







% --- Executes just before matlab_visualization is made visible.
function matlab_visualization_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to matlab_visualization (see VARARGIN)

% Choose default command line output for matlab_visualization
handles.output = hObject;


% Create a cell array (with only one element to start with) where to store
% all available information about different neighbours
handles.allNeighbours{1}.s = str2double(get(handles.edit_s,'String'));
handles.allNeighbours{1}.L = str2double(get(handles.edit_L,'String'));
handles.allNeighbours{1}.r = str2double(get(handles.edit_r,'String'));
handles.allNeighbours{1}.theta = get(handles.slider_theta,'Value')*pi/180;
handles.allNeighbours{1}.phi = get(handles.slider_phi,'Value')*pi/180;

% Initialize two attraction points around each neighbour
handles.allNeighbours{1}.d(1) = str2double(get(handles.edit_d, 'String'));
handles.allNeighbours{1}.alpha(1) = (mod(180 - str2double(get(handles.edit_alpha, 'String')), 360) - 180)*pi/180;
handles.allNeighbours{1}.d(2) = str2double(get(handles.edit_d, 'String'));
handles.allNeighbours{1}.alpha(2) = str2double(get(handles.edit_alpha, 'String'))*pi/180;

set(handles.popupmenu_selected_attractor, 'Max', 2);

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using matlab_visualization.
if strcmp(get(hObject,'Visible'),'off')
    launchUpdate(hObject, eventdata, handles);
end








% --- Outputs from this function are returned to the command line.
function varargout = matlab_visualization_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;







% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % --------------------------------------------------------------------
% function OpenMenuItem_Callback(hObject, eventdata, handles)
% % hObject    handle to OpenMenuItem (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% file = uigetfile('*.fig');
% if ~isequal(file, 0)
%     open(file);
% end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
    ['Close ' get(handles.figure1,'Name') '...'],...
    'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)




% --- Executes on slider movement.
function sliderPhi_Callback(hObject, eventdata, handles) % slider1 controls phi of the neighbour
% hObject    handle to slider_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% disp(sprintf('New value : %f', get(hObject, 'Value')));
set(handles.edit_phi, 'String', num2str(get(hObject,'Value')));

% assign the new phi value to the selected neighbour
currObj = get(handles.popupmenu_selected_neighbour, 'Value');
handles.allNeighbours{currObj}.phi = get(hObject,'Value')*pi/180;

% Update handles structure
guidata(hObject, handles);

launchUpdate(hObject, eventdata, handles);







% --- Executes during object creation, after setting all properties.
function sliderPhi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end









function edit_phi_Callback(hObject, eventdata, handles) % edit1 controls phi of the neighbour
% hObject    handle to edit_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_phi as text
%        str2double(get(hObject,'String')) returns contents of edit_phi as a double
newValue = str2double(get(hObject,'String'));
if (newValue > 180)
    set(hObject,'String', num2str(180));
end
if (newValue < -180)
    set(hObject,'String', num2str(-180));
end

set(handles.slider_phi, 'Value', str2double(get(hObject,'String')));
% assign the new phi value to the selected neighbour
currObj = get(handles.popupmenu_selected_neighbour, 'Value');
handles.allNeighbours{currObj}.phi = str2double(get(hObject,'String'))*pi/180;

% Update handles structure
guidata(hObject, handles);

launchUpdate(hObject, eventdata, handles);






% --- Executes during object creation, after setting all properties.
function edit_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end











% --- Executes on slider movement.
function slider_theta_Callback(hObject, eventdata, handles) % slider2 controls theta of the neighbour
% hObject    handle to slider_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.edit_theta, 'String', num2str(get(hObject,'Value')));

% assign the new phi value to the selected neighbour
currObj = get(handles.popupmenu_selected_neighbour, 'Value');
handles.allNeighbours{currObj}.theta = get(hObject,'Value')*pi/180;

% Update handles structure
guidata(hObject, handles);

launchUpdate(hObject, eventdata, handles);








% --- Executes during object creation, after setting all properties.
function slider_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end












function edit_theta_Callback(hObject, eventdata, handles) % edit_theta controls theta of the neighbour
% hObject    handle to edit_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_theta as text
%        str2double(get(hObject,'String')) returns contents of edit_theta as a double
newValue = str2double(get(hObject,'String'));
if (newValue > 180)
    set(hObject,'String', num2str(180));
end
if (newValue < -180)
    set(hObject,'String', num2str(-180));
end


set(handles.slider_theta, 'Value', str2double(get(hObject,'String')));

currObj = get(handles.popupmenu_selected_neighbour, 'Value');
handles.allNeighbours{currObj}.theta = str2double(get(hObject,'String'))*pi/180;

% Update handles structure
guidata(hObject, handles);

launchUpdate(hObject, eventdata, handles);










% --- Executes during object creation, after setting all properties.
function edit_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end













function edit_r_Callback(hObject, eventdata, handles) % edit_r controls "r" (the distance to the neighbour)
% hObject    handle to edit_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_r as text
%        str2double(get(hObject,'String')) returns contents of edit_r as a double

currObj = get(handles.popupmenu_selected_neighbour, 'Value');
handles.allNeighbours{currObj}.r = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

launchUpdate(hObject, eventdata, handles);













% --- Executes during object creation, after setting all properties.
function edit_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end















function edit_L_Callback(hObject, eventdata, handles) % edit_L controls "L", the length of the neighbour
% hObject    handle to edit_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_L as text
%        str2double(get(hObject,'String')) returns contents of edit_L as a double


currObj = get(handles.popupmenu_selected_neighbour, 'Value');
handles.allNeighbours{currObj}.L = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

launchUpdate(hObject, eventdata, handles);













% --- Executes during object creation, after setting all properties.
function edit_L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end









function edit_sf_Callback(hObject, eventdata, handles) % edit5 controls "sf", the speed of the focal fish
% hObject    handle to edit_sf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sf as text
%        str2double(get(hObject,'String')) returns contents of edit_sf as a double
launchUpdate(hObject, eventdata, handles);










% --- Executes during object creation, after setting all properties.
function edit_sf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end









function edit_s_Callback(hObject, eventdata, handles) % edit6 controls "s", the speed of the neighbour
% hObject    handle to edit_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_s as text
%        str2double(get(hObject,'String')) returns contents of edit_s as a double


currObj = get(handles.popupmenu_selected_neighbour, 'Value');
handles.allNeighbours{currObj}.s = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

launchUpdate(hObject, eventdata, handles);








% --- Executes during object creation, after setting all properties.
function edit_s_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end








function edit_d_Callback(hObject, eventdata, handles) % edit7 controls the intensity of the background flow
% hObject    handle to edit_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_d as text
%        str2double(get(hObject,'String')) returns contents of edit_d as a double


currAttr = get(handles.popupmenu_selected_attractor, 'Value');
handles.allNeighbours{:}.d(currAttr) = str2double(get(hObject,'String'));


% Update handles structure
guidata(hObject, handles);

launchUpdate(hObject, eventdata, handles);








% --- Executes during object creation, after setting all properties.
function edit_d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







function edit8_Callback(hObject, eventdata, handles) % edit8 controls the spatial frequency of the filter
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
launchUpdate(hObject, eventdata, handles);




% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles) % edit9 controls the width of the frequency filter
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
launchUpdate(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on selection change in popupmenu_selected_neighbour.
function popupmenu_selected_neighbour_Callback(hObject, eventdata, handles) % popupmenu2 changes the current neighbour
% hObject    handle to popupmenu_selected_neighbour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_selected_neighbour contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_selected_neighbour
currObj = get(hObject, 'Value');

set(handles.edit_s,'String', num2str(handles.allNeighbours{currObj}.s));
set(handles.edit_L,'String', num2str(handles.allNeighbours{currObj}.L));
set(handles.edit_r,'String', num2str(handles.allNeighbours{currObj}.r));
set(handles.edit_theta,'String', num2str(handles.allNeighbours{currObj}.theta*180/pi));
set(handles.edit_phi,'String', num2str(handles.allNeighbours{currObj}.phi*180/pi));
set(handles.slider_theta, 'Value', (handles.allNeighbours{currObj}.theta)*180/pi);
set(handles.slider_phi,'Value', (handles.allNeighbours{currObj}.phi)*180/pi);









% --- Executes during object creation, after setting all properties.
function popupmenu_selected_neighbour_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_selected_neighbour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end








% --- Executes on button press in pushbutton_add_neighbour.
function pushbutton_add_neighbour_Callback(hObject, eventdata, handles) % pushbuttonAddNeighbour adds a new neighbour
% hObject    handle to pushbutton_add_neighbour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currMaxNObjects = get(handles.popupmenu_selected_neighbour, 'Max');
currMaxNObjects = currMaxNObjects + 1;
set(handles.popupmenu_selected_neighbour, 'Max', currMaxNObjects);
set(handles.popupmenu_selected_neighbour, 'String', 1:currMaxNObjects);


currObj = currMaxNObjects;
set(handles.popupmenu_selected_neighbour, 'Value', currObj);

handles.allNeighbours{currObj}.s = 1.0;
handles.allNeighbours{currObj}.L = 1.0;
handles.allNeighbours{currObj}.r = 2+rand(1)*3;
handles.allNeighbours{currObj}.theta = (rand(1)*360-180)*pi/180;
handles.allNeighbours{currObj}.phi = (rand(1)*90-45)*pi/180;
handles.allNeighbours{currObj}.d = handles.allNeighbours{currObj - 1}.d;
handles.allNeighbours{currObj}.alpha = handles.allNeighbours{currObj - 1}.alpha;

set(handles.edit_s,'String', num2str(handles.allNeighbours{currObj}.s));
set(handles.edit_L,'String', num2str(handles.allNeighbours{currObj}.L));
set(handles.edit_r,'String', num2str(handles.allNeighbours{currObj}.r));
set(handles.edit_theta,'String', num2str(handles.allNeighbours{currObj}.theta*180/pi));
set(handles.edit_phi,'String', num2str(handles.allNeighbours{currObj}.phi*180/pi));
set(handles.slider_theta, 'Value', (handles.allNeighbours{currObj}.theta)*180/pi);
set(handles.slider_phi,'Value', (handles.allNeighbours{currObj}.phi)*180/pi);


% Update handles structure
guidata(hObject, handles);

launchUpdate(hObject, eventdata, handles);



% --- Executes on mouse press over axes background.
function drawing_ButtonDownFcn(hObject, eventdata, handles) % this has no effect; I am only interested in button down responses to objects, not to the background
% hObject    handle to drawing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pt = get(handles.drawing, 'CurrentPoint');
disp('click at point:');
pt
% get(gcbo,'Tag')











function launchUpdate(hObject, eventdata, handles)


sff = str2double(get(handles.edit_sf,'String')); % speed of the focal fish

nNeigh = length(handles.allNeighbours); % number of neighbours
nAttr = length(handles.allNeighbours{1}.d); % the number of attractors is the same for all neighbours

% for convenience I make arrays for each parameter of neighbour's body and
% position
for jj = 1:nNeigh % loop on neighbours
    s(jj) = handles.allNeighbours{jj}.s; % speed of neighbour
    L(jj) = handles.allNeighbours{jj}.L; % length of neighbour
    r(jj) = handles.allNeighbours{jj}.r; % distance of neighbour
    theta(jj) = handles.allNeighbours{jj}.theta; % bearing of neighbour in the frame of reference of focal
    phi(jj) = handles.allNeighbours{jj}.phi; % difference of bearing between neighbour and focal
    xNeigh(jj) = r(jj)*cos(theta(jj)); % x coordinate of the neighbour, relative to focal
    yNeigh(jj) = r(jj)*sin(theta(jj)); % y coordinate of neighbour, relative to focal
    for ii = 1:length(handles.allNeighbours{jj}.d) % loop on attractors
        d(jj,ii) = handles.allNeighbours{jj}.d(ii); % distance of the attractor from the body of the neighbour
        alpha(jj,ii) = handles.allNeighbours{jj}.alpha(ii); % angle of the attractor relative to the moving direction of the neighbour
        xAttr(jj,ii) = xNeigh(jj) + d(jj,ii)*cos(alpha(jj,ii) + phi(jj)); % x position of the attractor
        yAttr(jj,ii) = yNeigh(jj) + d(jj,ii)*sin(alpha(jj,ii) + phi(jj)); % y position of the attractor
    end
end


%% In handles.drawing draws the neighbours and their attractors
axes(handles.drawing);
cla;
hold on;
% line([-L, 0], [0, 0], 'LineWidth', 2);
plot(0,0, 'Marker', 'o', 'MarkerSize', 12);
focalSpeedHandle = quiver(0,0,0,sff,0);

myColourMap = jet(max(1,length(handles.allNeighbours)));


for jj = 1: nNeigh
    hN = line([r(jj)*sin(theta(jj)) - L(jj)/2*sin(phi(jj)), r(jj)*sin(theta(jj)) + L(jj)/2*sin(phi(jj))], ...
        [r(jj)*cos(theta(jj)) - L(jj)/2*cos(phi(jj)), r(jj)*cos(theta(jj)) + L(jj)/2*cos(phi(jj))], ...
        'LineWidth', 2, 'Color', myColourMap(jj,:));
    set(hN,'ButtonDownFcn',{@mouseClickCallback,handles}, 'Tag', num2str(jj));
    % set(handleNeighbourBody, 'Tag', num2str(get(handles.popupmenu_selected_neighbour, 'Value')));
    neighbourSpeedHandle = quiver(r(jj)*sin(theta(jj)) + L(jj)/2*sin(phi(jj)), r(jj)*cos(theta(jj)) + L(jj)/2*cos(phi(jj)), ...
        s(jj)*sin(phi(jj)), s(jj)*cos(phi(jj)), 0);
    set(neighbourSpeedHandle, 'Color', myColourMap(jj,:));
    
    for ii = 1:nAttr
        % draw the current position of the attractor
        plot(yAttr(jj,ii), xAttr(jj,ii), 'Marker', 'x', 'MarkerSize', 14, 'Color', 'k', 'LineStyle', 'none');
        
        % draw the position of the attractor at the matching time
        matchingTime = 1;
        plot(yAttr(jj,ii)  + [0,matchingTime] * s(jj)*sin(phi(jj)), xAttr(jj,ii) + [0,matchingTime] * s(jj)*cos(phi(jj)), 'Marker', 'none', 'MarkerSize', 14, 'Color', 'k', 'LineStyle', '-');
    end
    
    
end


% find the closest attractor (the one to which the focal is actually attracted)
[~, closestAttr] = min(xAttr(:).^2 + yAttr(:).^2);
[jj,ii] = ind2sub([nNeigh, nAttr], closestAttr);

% mark the closest attractor in red
plot(yAttr(jj,ii), xAttr(jj,ii), 'Marker', 'x', 'MarkerSize', 14, 'Color', 'r', 'LineStyle', 'none');


% draw the position of the closest attractor at the matching time
matchingTime = 1;
plot(yAttr(jj,ii)  + [0,matchingTime] * s(jj)*sin(phi(jj)), xAttr(jj,ii) + [0, matchingTime] * s(jj)*cos(phi(jj)), 'Marker', 'none', 'MarkerSize', 14, 'Color', 'r', 'LineStyle', '-');

plot([0,yAttr(closestAttr) + matchingTime * s(jj)*sin(phi(jj))], [sff, xAttr(closestAttr) + matchingTime* s(jj)*cos(phi(jj))], 'Color', 'b', 'LineStyle', ':', 'LineWidth', 1.2) % plot a line between the focal individual and its attractor




function mouseClickCallback(hObject, ~, handles)
currObj = str2double(get(hObject, 'Tag'));
set(handles.popupmenu_selected_neighbour, 'Value', currObj);

set(handles.edit_s,'String', num2str(handles.allNeighbours{currObj}.s));
set(handles.edit_L,'String', num2str(handles.allNeighbours{currObj}.L));
set(handles.edit_r,'String', num2str(handles.allNeighbours{currObj}.r));
set(handles.edit_theta,'String', num2str(handles.allNeighbours{currObj}.theta*180/pi));
set(handles.edit_phi,'String', num2str(handles.allNeighbours{currObj}.phi*180/pi));
set(handles.slider_theta, 'Value', (handles.allNeighbours{currObj}.theta)*180/pi);
set(handles.slider_phi,'Value', (handles.allNeighbours{currObj}.phi)*180/pi);


% retinalSize = sum(diff(thetaI));
% translation = mean(opticFlow) * abs(retinalSize);
% expansion = sum(diff(opticFlow)) *sign(retinalSize);
% % set(handles.text11, 'String', num2str(translation));
% % set(handles.text10, 'String', num2str(expansion));
% % set(handles.text13, 'String', num2str(retinalSize));
% %









function fourierFilter = draw_1d_fft_filter(x0, sigma, len) % create the fourier filter to analyse the optic flow
X = 1:len;
fourierFilter = ifftshift(exp(-(log(abs((X-round(len/2-1))./(x0))).^2)/(2*sigma*sigma)));
fourierFilter(1)=0; % this puts the DC to zero
fourierFilter(round(len/2+1):len) = 0; % this allows to compute the even symmetric and odd symmetric components in one single step
% figure, plot(imag(fftshift(ifft(fourierFilter)))) % odd symmetric
% figure, plot(real(fftshift(ifft(fourierFilter)))) % even symmetric



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_acceleration_Callback(hObject, eventdata, handles)
% hObject    handle to edit_acceleration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_acceleration as text
%        str2double(get(hObject,'String')) returns contents of edit_acceleration as a double


% --- Executes during object creation, after setting all properties.
function edit_acceleration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_acceleration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_add_attractor.
function pushbutton_add_attractor_Callback(hObject, eventdata, handles) % this adds a new attractor
% hObject    handle to pushbutton_add_attractor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currMaxNAttractors = get(handles.popupmenu_selected_attractor, 'Max');
currMaxNAttractors = currMaxNAttractors + 1;
set(handles.popupmenu_selected_attractor, 'Max', currMaxNAttractors);
set(handles.popupmenu_selected_attractor, 'String', 1:currMaxNAttractors);


currAttr = currMaxNAttractors;
set(handles.popupmenu_selected_attractor, 'Value', currAttr);

handles.allNeighbours{:}.d(currAttr) = 2+rand(1)*3;
handles.allNeighbours{:}.alpha(currAttr) = (rand(1)*360-180)*pi/180;

set(handles.edit_d,'String', num2str(handles.allNeighbours{1}.d(currAttr)));
set(handles.edit_alpha,'String', num2str(handles.allNeighbours{1}.alpha(currAttr)*180/pi));



% Update handles structure
guidata(hObject, handles);

launchUpdate(hObject, eventdata, handles);





function edit_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newValue = str2double(get(hObject,'String'));
if (newValue > 180)
    set(hObject,'String', num2str(180));
end
if (newValue < -180)
    set(hObject,'String', num2str(-180));
end

currAttr = get(handles.popupmenu_selected_attractor, 'Value');
handles.allNeighbours{:}.alpha(currAttr) = str2double(get(hObject,'String'))*pi/180;

% Update handles structure
guidata(hObject, handles);

launchUpdate(hObject, eventdata, handles);



% Hints: get(hObject,'String') returns contents of edit_alpha as text
%        str2double(get(hObject,'String')) returns contents of edit_alpha as a double


% --- Executes during object creation, after setting all properties.
function edit_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_selected_attractor.
function popupmenu_selected_attractor_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_selected_attractor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_selected_attractor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_selected_attractor

currAttr = get(hObject, 'Value');

set(handles.edit_d,'String', num2str(handles.allNeighbours{1}.d(currAttr)));
set(handles.edit_alpha,'String', num2str(handles.allNeighbours{1}.alpha(currAttr)*180/pi));






% --- Executes during object creation, after setting all properties.
function popupmenu_selected_attractor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_selected_attractor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_update_plots.
function pushbutton_update_plots_Callback(hObject, eventdata, handles)
% function pushbutton_update_plots_Callback
% this function updates the four heatmap plots on the right of the figure
% the plots express turning and acceleration responses as a function of
% parameters
% Plot top-left: turning response as a function of r (distance to the
% neighbour) and theta (bearing to the neighbour).


% get selected neighbour
matchingTime = 1;
currObj = get(handles.popupmenu_selected_neighbour, 'Value');
sff = str2double(get(handles.edit_sf,'String')); % speed of the focal fish
s = handles.allNeighbours{currObj}.s;
nAttr = length(handles.allNeighbours{currObj}.d); % number of attractors
[useFontSize, myColourMap] = set_figure_parameters();




% bins for the figures
rBins = 0:0.5:10;
% thetaBins = linspace(-pi + pi/3, pi - pi/3, 25); % radians (13 bins make 30 degrees each)
thetaBins = linspace(-pi, pi, 25); % radians (13 bins make 30 degrees each)

centerThetaBins = thetaBins(1:end-1) + diff(thetaBins)/2;
centerRBins = rBins(1:end-1) + diff(rBins)/2;

useAllPhiBins = 0;
if useAllPhiBins
    phiBins = linspace(-pi, pi, 25); % radians (13 bins make 30 degrees each)
    centerPhiBins = phiBins(1:end-1) + diff(phiBins)/2;
else
    phiBins = linspace(-pi/20, pi/20, 11);
    centerPhiBins = phiBins(1:end-1) + diff(phiBins)/2;
end


% initialize empty arrays
turningRTheta = NaN(length(rBins), length(thetaBins));
accelerationRTheta = NaN(length(rBins), length(thetaBins));

turningPhiTheta = NaN(length(phiBins), length(thetaBins));
accelerationPhiTheta = NaN(length(phiBins), length(thetaBins));

defaultPhi = handles.allNeighbours{currObj}.phi;
defaultR = handles.allNeighbours{currObj}.r;


for countTheta=1:length(thetaBins)-1
    currTheta = centerThetaBins(countTheta);
    
    
    for countR=1:length(rBins)-1
        currR = centerRBins(countR);
        
        
        % closest attractor
        d = nan(nAttr,1);
        alpha = nan(nAttr,1);
        xAttr = nan(nAttr,1);
        yAttr = nan(nAttr,1);
        for ii = 1:nAttr % loop on attractors
            d(ii) = handles.allNeighbours{currObj}.d(ii); % distance of the attractor from the body of the neighbour
            alpha(ii) = handles.allNeighbours{currObj}.alpha(ii); % angle of the attractor relative to the moving direction of the neighbour
            xAttr(ii) = currR*cos(currTheta) + d(ii)*cos(alpha(ii) + defaultPhi); % x position of the attractor at time 0
            yAttr(ii) = currR*sin(currTheta) + d(ii)*sin(alpha(ii) + defaultPhi); % y position of the attractor at time 0
        end
        
        [~, closestAttr] = min(xAttr.^2 + yAttr.^2);
        
        turningRTheta(countR, countTheta) = sign(yAttr(closestAttr) + matchingTime * s * sin(defaultPhi));
        accelerationRTheta(countR, countTheta) = sign(xAttr(closestAttr) + matchingTime * s * cos(defaultPhi) -sff);
    end
    
    for countPhi=1:length(phiBins)-1
        
        currPhi = centerPhiBins(countPhi);
        
        % closest attractor
        d = nan(nAttr,1);
        alpha = nan(nAttr,1);
        xAttr = nan(nAttr,1);
        yAttr = nan(nAttr,1);
        for ii = 1:nAttr % loop on attractors
            d(ii) = handles.allNeighbours{currObj}.d(ii); % distance of the attractor from the body of the neighbour
            alpha(ii) = handles.allNeighbours{currObj}.alpha(ii); % angle of the attractor relative to the moving direction of the neighbour
            xAttr(ii) = defaultR*cos(currTheta) + d(ii)*cos(alpha(ii) + currPhi); % x position of the attractor at time 0
            yAttr(ii) = defaultR*sin(currTheta) + d(ii)*sin(alpha(ii) + currPhi); % y position of the attractor at time 0
        end
        
        [~, closestAttr] = min(xAttr.^2 + yAttr.^2);
        
        turningPhiTheta(countPhi, countTheta) = sign(yAttr(closestAttr) + matchingTime * s * sin(currPhi));
        accelerationPhiTheta(countPhi, countTheta) = sign(xAttr(closestAttr) + matchingTime * s * cos(currPhi) -sff);
    end
end


X = rBins'*cos(pi/2 -thetaBins); % I add pi/2 because I want the polar axes with theta = 0 on top; I subtract thetaBins because I want clockwise axes
Y = rBins'*sin(pi/2 -thetaBins);









%% figure for turning angle
valueForColourScale = 0; % set a maximum turning angle at 30 degrees
axes(handles.turning_r_theta);
cla;
%    title('turning vs. r and theta');
pcolor(X,Y,turningRTheta*180/pi); caxis([-valueForColourScale, valueForColourScale]);
hold on;
fill([0, -0.5, 0, 0.5],[1.0, -0.8, 0, -0.8], 'r', 'LineWidth', 1, 'LineStyle', 'none');
axis equal xy tight;
box off;
colormap(myColourMap);

set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(rBins):1:max(rBins), 'YTick', -max(rBins):1:max(rBins), ...
    'XLim', [-max(rBins)*12/10 max(rBins)*12/10], 'YLim', [-max(rBins)*12/10, max(rBins)*12/10], ... 'XTickLabel', -max(rBins):1:max(rBins), 'YTickLabel',  -max(rBins):1:max(rBins), ...
    'FontName', 'Arial');
text(0,max(rBins)*11/10,'\vartheta = 0', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(max(rBins)*11/10,0,'\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(-max(rBins)*11/10,0,'-\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(0,-max(rBins)*11/10,'\pi', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');

plot((d).*sin((-alpha + defaultPhi)) - s*sin(defaultPhi), (d).*cos(-alpha + defaultPhi) - s*cos(defaultPhi) + sff, 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2, 'LineStyle', 'none', 'MarkerFaceColor', [0.5 0.5 0.5], 'Color', [0 0 0]);


%     % draw the position where acceleration and turning are zero
%     matchingTime = 1;
%     xAttr(jj,ii) = xNeigh(jj) + d(jj,ii)*cos(alpha(jj,ii) + phi(jj)); % x position of the attractor
%     yAttr(jj,ii) = yNeigh(jj) + d(jj,ii)*sin(alpha(jj,ii) + phi(jj)); % y position of the attractor
%
%     plot(yAttr(jj,ii)  + [0,matchingTime] * s(jj)*sin(phi(jj)), xAttr(jj,ii) + [0, matchingTime] * s(jj)*cos(phi(jj)), 'Marker', 'none', 'MarkerSize', 14, 'Color', 'r', 'LineStyle', '-');




xlabel('distance from neighbour'); ylabel('distance from neighbour');
title(sprintf('Neighbour orientation phi=%0.1f', defaultPhi), 'FontName', 'Arial');
%     colorbarHandle = colorbar;
%     set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');
%     ylabel(colorbarHandle, 'turning angle', 'FontSize', useFontSize);
colormap(myColourMap);

xLim = get(gca, 'XLim');
yLim = get(gca, 'YLim');

text((xLim(2) - xLim(1))*0.1 + xLim(1), (yLim(2) - yLim(1))*0.9 + yLim(1),'turning', 'FontWeight', 'bold');






%% figure for acceleration
valueForColourScale = 0;
axes(handles.acceleration_r_theta);
cla;
%    title('acceleration vs. r and theta');
pcolor(X,Y,accelerationRTheta); caxis([-valueForColourScale, valueForColourScale]);
hold on;
fill([0, -0.5, 0, 0.5],[1.0, -0.8, 0, -0.8], 'r', 'LineWidth', 1, 'LineStyle', 'none');
axis equal xy tight;
box off;
colormap(myColourMap);

set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(rBins):1:max(rBins), 'YTick', -max(rBins):1:max(rBins), ...
    'XLim', [-max(rBins)*12/10 max(rBins)*12/10], 'YLim', [-max(rBins)*12/10, max(rBins)*12/10], ... 'XTickLabel', -max(rBins):1:max(rBins), 'YTickLabel',  -max(rBins):1:max(rBins), ...
    'FontName', 'Arial');
text(0,max(rBins)*11/10,'\vartheta = 0', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(max(rBins)*11/10,0,'\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(-max(rBins)*11/10,0,'-\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(0,-max(rBins)*11/10,'\pi', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');

% plot(d.*cos(pi*3/2 - alpha - defaultPhi), d.*sin(pi*3/2 - alpha - defaultPhi), 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2, 'LineStyle', 'none', 'MarkerFaceColor', [0.5 0.5 0.5], 'Color', [0 0 0]);
plot((d).*sin((-alpha + defaultPhi)) - s*sin(defaultPhi), (d).*cos(-alpha + defaultPhi) - s*cos(defaultPhi) + sff, 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2, 'LineStyle', 'none', 'MarkerFaceColor', [0.5 0.5 0.5], 'Color', [0 0 0]);

% text(-max(rBins)*11/10, max(rBins)*11/10, 'C', 'FontName', 'Arial', 'FontSize', 50);

%
% if defaultPhi == 0
% text(-max(rBins)*11/10, +max(rBins)*11/10, {'\phi=0\pi'}, 'FontName', 'Arial', 'FontSize', useFontSize, 'Color', 'k');
% else
% text(-max(rBins)*11/10, +max(rBins)*11/10, {'\phi =', num2str(rats(defaultPhi/pi)), '\pi'}, 'FontName', 'Arial', 'FontSize', useFontSize, 'Color', 'k');
% end


xlabel('distance from neighbour'); % ylabel('distance from neighbour');
title(sprintf('Neighbour orientation phi=%0.1f', defaultPhi), 'FontName', 'Arial');
%     colorbarHandle = colorbar;
%     set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
%     ylabel(colorbarHandle, 'acceleration', 'FontSize', useFontSize);
colormap(myColourMap);


xLim = get(gca, 'XLim');
yLim = get(gca, 'YLim');

text((xLim(2) - xLim(1))*0.1 + xLim(1), (yLim(2) - yLim(1))*0.9 + yLim(1),'acceleration', 'FontWeight', 'bold');








% figure for turning angle phi theta
valueForColourScale = (max(max(abs(turningPhiTheta))));
axes(handles.turning_phi_theta);
cla;
%    title('turning vs. phi and theta');
imagesc(turningPhiTheta(1:end-1, 1:end-1)*180/pi , [-valueForColourScale, valueForColourScale]);


if useAllPhiBins
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:4:length(thetaBins), 'YTick', ...
        0.5:4:length(phiBins), 'TickDir', 'out', 'XTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
        'YTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
        'FontName', 'symbol');
else
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:4:length(thetaBins), 'YTick', ...
        linspace(0.5, length(phiBins)-0.5, 7), 'TickDir', 'out', 'XTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
        'YTickLabel', {'-p/20', '-p/30', '-p/60', '0', 'p/60', 'p/30', 'p/20'}, 'FontName', 'symbol');
end


xlabel('\theta (rad)', 'FontName', 'arial'); ylabel('\phi (rad)', 'FontName', 'arial');
%     colorbarHandle = colorbar;
%     set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
%     ylabel(colorbarHandle, 'turning (rad / time step)', 'FontSize', useFontSize);
axis xy equal tight;

title(sprintf('Neighbour distance r=%0.1f', defaultR), 'FontName', 'Arial');
colormap(myColourMap);


xLim = get(gca, 'XLim');
yLim = get(gca, 'YLim');

text((xLim(2) - xLim(1))*0.1 + xLim(1), (yLim(2) - yLim(1))*0.9 + yLim(1),'turning', 'FontWeight', 'bold', 'Color', 'r');








% figure for acceleration phi theta
valueForColourScale = (max(max(abs(accelerationPhiTheta))));
axes(handles.acceleration_phi_theta);
cla;
%    title('acceleration vs. phi and theta');

imagesc(accelerationPhiTheta(1:end-1, 1:end-1) , [-valueForColourScale, valueForColourScale]);

if useAllPhiBins
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:4:length(thetaBins), 'YTick', ...
        0.5:4:length(phiBins), 'TickDir', 'out', 'XTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
        'YTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
        'FontName', 'symbol');
else
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:4:length(thetaBins), 'YTick', ...
        linspace(0.5, length(phiBins)-0.5, 7), 'TickDir', 'out', 'XTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
        'YTickLabel', {'-p/20', '-p/30', '-p/60', '0', 'p/60', 'p/30', 'p/20'}, 'FontName', 'symbol');
end

xlabel('\theta (rad)', 'FontName', 'arial'); ylabel('\phi (rad)', 'FontName', 'arial');
%     colorbarHandle = colorbar;
%     set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
%     ylabel(colorbarHandle, 'acceleration (cm/s^2)', 'FontSize', useFontSize);
axis xy equal tight;
title(sprintf('Neighbour distance r=%0.1f', defaultR), 'FontName', 'Arial');
colormap(myColourMap);


xLim = get(gca, 'XLim');
yLim = get(gca, 'YLim');

text((xLim(2) - xLim(1))*0.1 + xLim(1), (yLim(2) - yLim(1))*0.9 + yLim(1),'acceleration', 'FontWeight', 'bold', 'Color', 'r');


launchUpdate(hObject, eventdata, handles);

[stablePointsTheta, stablePointsRho] = cart2pol((d).*cos((- alpha + defaultPhi)) - s*cos(defaultPhi)+sff, (d).*sin(- alpha + defaultPhi) - s*sin(defaultPhi));
stablePointsTheta = (mod(pi + stablePointsTheta, 2*pi) - pi) * 180/pi;

axes(handles.drawing);
xLim = get(gca, 'XLim');
yLim = get(gca, 'YLim');

stablePointsString = {'Stable points:'};
for jj = 1:length(stablePointsTheta)
    newString = sprintf('theta = %f; rho = %f', stablePointsTheta(jj), stablePointsRho(jj));
    stablePointsString{end+1} = newString;
end
text((xLim(2) - xLim(1))*0.1 + xLim(1), (yLim(2) - yLim(1))*0.9 + yLim(1),stablePointsString, 'FontWeight', 'bold');







function [useFontSize, myColourMap] = set_figure_parameters()
% Function set_figure_parameters
% this function sets some figure parameters, such as the font size for
% figure annotations and the colourmap

useFontSize = 12; % font size for the figures

% create a colour blind friendly colormap
x=(0:1:300)/300;

R = x/0.32 - 0.78125;
R(R<0) = 0;
R(R>1) = 1;

G = 2*x - 0.84;
G(G<0) = 0;
G(G>1) = 1;


B(find(x<=0.25)) = 4*x(x<=0.25);
B(intersect(find(x>0.25), find(x<=0.42))) = 1;
B(intersect(find(x>0.42), find(x<=0.92))) = -2*x(x>0.42 & x<=0.92) + 1.84;
B(find(x>0.92)) = x(x>0.92)/0.08 - 11.5;

myColourMap=[R', G', B'];

myColourMap=myColourMap(1:256,:); % I want to end with yellow

% here use instead a standard colormap
% myColourMap = cool(256);% this line actua



% --- Executes on button press in pushbutton_rotate_90_degrees.
function pushbutton_rotate_90_degrees_Callback(hObject, eventdata, handles)
% Function rotate_90_degrees_Callback
% this function reads the angular positions describing the attractors
% and rotates them by + 90 degrees


nAttr = length(handles.allNeighbours{1}.d); % number of attractors (the number is the same for all neighbours)
currAttr = get(handles.popupmenu_selected_attractor, 'Value'); % current (selected) attractor

for jj = 1:nAttr
    handles.allNeighbours{:}.alpha(jj) = (mod(pi + handles.allNeighbours{1}.alpha(jj) + pi/2, 2*pi) - pi); % 90 degrees rotation
end

set(handles.edit_alpha,'String', num2str(handles.allNeighbours{1}.alpha(currAttr)*180/pi));
% Update handles structure
guidata(hObject, handles);

launchUpdate(hObject, eventdata, handles);
