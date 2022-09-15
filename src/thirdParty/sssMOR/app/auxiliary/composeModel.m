function varargout = composeModel(varargin)
% GUI to create state space models from (1st or 2nd order) system matrices
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Sylvia Cremer
% Last Change:  13 Feb 2011
% ------------------------------------------------------------------



%Structure:
%
%   GUI INITIALIZATION
%
%   FIGURE(TOP-LEVEL)
%
%       TOP MENUE
%
%       PANEL ABCD
%
%       PANEL MDK
%
%       BOTTOM MENUE
%
%   UTILITY FUNCTIONS






%--------------------------------------------------------------------------
%                           GUI INITIALIZATION
%--------------------------------------------------------------------------

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @composeModel_OpeningFcn, ...
                   'gui_OutputFcn',  @composeModel_OutputFcn, ...
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

function composeModel_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
    % figure has been opened. evaluate parameters
    handles.output = hObject;
    handles.inputGUI=varargin{2};
    guidata(hObject, handles);
    set(handles.panel_ABCD,'Visible','on')
    set(handles.panel_MDK,'Visible','off')
    set(handles.pb_loadfile,'Visible','on')
    set(handles.pb_create,'Visible','on')
    set(handles.st_a,'Visible','off')
    set(handles.st_b,'Visible','off')
    set(handles.ed_a,'UserData',0)
    set(handles.ed_a,'UserData',0)
    list_matrices(handles)

    %Make latex code possible for all static text-fields whose tag names start
    %with "latex"  
    
    lbls = findall(hObject);
    
    for i=1:length(lbls)
        
          l = lbls(i);
          t = get(l,'tag');
          
          if strfind(t,'latex_') == 1
              
             % Get current text, position and tag
             t = 'nothing'; 
             set(l,'units','characters');
             s = get(l,'string');
             p = get(l,'position');
             parent = get(l,'parent');
             
             % Remove the UICONTROL
             delete(l);
             
             % Replace it with an axis with an text object on it
             handles.(t) = axes('parent',parent,'units','characters','position',[p(1) p(2) p(3) p(4)],'visible','off');
             handles.(t) = text(0,0.6,s,'interpreter','latex','FontSize',13);            
          end
    end
    
    %Preselect the option include matrices from workspace
    
    set(handles.cb_include,'Value',1);
    cb_include_Callback(hObject,eventdata,handles);
    
function varargout = composeModel_OutputFcn(hObject, eventdata, handles) %#ok<*INUSD>
varargout = {};



%--------------------------------------------------------------------------
%                            FIGURE(TOP LEVEL)
%--------------------------------------------------------------------------

function figure1_DeleteFcn(hObject, eventdata, handles)
    try
        set(handles.inputGUI.lb_systems,'Value',[])
        set(handles.inputGUI.lb_systems,'String',systemsInWorkspace)
        set(handles.inputGUI.lb_matrixes,'Value',[])
        set(handles.inputGUI.lb_matrixes,'String',matricesInWorkspace)
    catch exception
        if strcmp(exception.identifier,'MATLAB:class:InvalidHandle')
            % GUI has been closes, handles are not valid anymore
        else
            throw(exception)
        end
    end
    delete(handles.figure1)
    % clear variable that have been written to workspace
    evalin('base','clear load_var_dont_destroy')


    
%--------------------------------------------------------------------------
%                               TOP MENUE
%--------------------------------------------------------------------------

function pb_loadfile_Callback(hObject, eventdata, handles)
set(hObject,'Enable','off') % doppeltes aufrufen der uigetfile funktion vermeiden
temp1=sprintf('%s.mat',handles.inputGUI.letzterpfad);
[FileName,PathName] = uigetfile(temp1);
set(hObject,'Enable','on')
if FileName==0
    return
end
handles.inputGUI.letzterpfad=PathName;
guidata(hObject, handles);
if isempty(strfind(FileName,'.mat')) % nur laden von .mat Dateien m�glich
    errordlg('Only .mat files allowed','Error Dialog','modal')
    uiwait
    return
end 
% die in der Datei vorhandenen Daten werden zun�chst alle in den workspace
% geladen, um mehrmaligen Netzwerkzugriff zu vermeiden (Zeit)
evalin('base',sprintf('load_var_dont_destroy=load(''%s%s'');', PathName,FileName));
set(handles.cb_file,'Value',1)
list_matrices(handles)
set(handles.cb_file,'Enable','on')
s=sprintf('Selected file: %s',FileName);
set(handles.st_file,'String',s)
suggestion=suggestVarname(strcat('sys_',FileName(1:strfind(FileName,'.')-1)),[]);
if ~isempty(suggestion)
    set(handles.ed_sysname,'String',suggestion);
    set(handles.ed_sysname,'UserData',0)
end

function bg_ABCMDK_SelectionChangeFcn(hObject, eventdata, handles)
if eventdata.NewValue==handles.rb_ABCD
    % ABCD selected
    set(handles.panel_ABCD,'Visible','on')
    set(handles.panel_MDK,'Visible','off')
elseif eventdata.NewValue==handles.rb_MDK
    % MDK selected
    set(handles.panel_MDK,'Visible','on')
    set(handles.panel_ABCD,'Visible','off')
end

function cb_include_Callback(hObject, eventdata, handles)
% Matrizen aus dem workspace (nicht) mit auflisten, entweder mit Matrizen
% aus der geladenen Datei, oder ohne
list_matrices(handles)

function cb_file_Callback(hObject, eventdata, handles)
list_matrices(handles)



%--------------------------------------------------------------------------
%                             PANEL ABCD
%--------------------------------------------------------------------------

function lb_A_Callback(hObject, eventdata, handles)
%Deletes options from the lists for the other system matrices that are not
%consistent with the choise for A

%     %Get the selected matrix A
%     
%     A = get_matrix_from_ws(hObject);
%     
%     %Remove the inconsistent entries for all other matrices
%     
%     clearInconsistentEntries(handles.lb_B,handles.lists.B,size(A,1),[]);
%     clearInconsistentEntries(handles.lb_E,handles.lists.E,size(A,1),size(A,2));
%     clearInconsistentEntries(handles.lb_C,handles.lists.C,[],size(A,1));
    


%--------------------------------------------------------------------------
%                              PANEL MDK
%--------------------------------------------------------------------------

function ed_a_Callback(hObject, eventdata, handles)
    % check alpha
    
    try
        h=str2num(get(hObject,'String')); %#ok<ST2NM>
        if isempty(h)
            errordlg('Use ''.'' as decimal seperator. Don''t use characters','Error Dialog','modal')
            uiwait
            set(hObject,'UserData',1)
        elseif imag(h)~=0
            errordlg('No imaginary numbers allowed','Error Dialog','modal')
            uiwait
            set(hObject,'UserData',1)
        elseif length(h)>1
            if length(h)==2
                h=sprintf('%i.%i',round(h(1)),round(h(2)));
                set(hObject,'String',h)
                set(hObject,'UserData',0)
                return
            end
            errordlg('Use ''.'' as decimal seperator','Error Dialog','modal')
            uiwait
            set(hObject,'UserData',1)
        else
            set(hObject,'UserData','0')
        end
    end

function ed_b_Callback(hObject, eventdata, handles)
    % check beta
    
    try
        h=str2num(get(hObject,'String')); %#ok<ST2NM>
        if isempty(h)
            errordlg('Use ''.'' as decimal seperator. Don''t use characters','Error Dialog','modal')
            uiwait
            set(hObject,'UserData',1)
        elseif imag(h)~=0
            errordlg('No imaginary numbers allowed','Error Dialog','modal')
            uiwait
            set(hObject,'UserData',1)
        elseif length(h)>1
            if length(h)==2
                h=sprintf('%i.%i',round(h(1)),round(h(2)));
                set(hObject,'String',h)
                set(hObject,'UserData',0)
                return
            end
            errordlg('Use ''.'' as decimal seperator','Error Dialog','modal')
            uiwait
            set(hObject,'UserData',1)
        else
            set(hObject,'UserData','0')
        end
    end
    
function lb_Mmdk_Callback(hObject, eventdata, handles)
%Deletes options from the lists for the other system matrices that are not
%consistent with the choise for M

    %Get the selected matrix M
    
%     M = get_matrix_from_ws(hObject);
%     
%     %Remove the inconsistent entries for all other matrices
%     
%     clearInconsistentEntries(handles.lb_Dmdk,handles.lists.Dmdk,size(M,1),size(M,1));
%     clearInconsistentEntries(handles.lb_Kmdk,handles.lists.Kmdk,size(M,1),size(M,1));
%     clearInconsistentEntries(handles.lb_Bmdk,handles.lists.Bmdk,size(M,1),[]);
%     clearInconsistentEntries(handles.lb_Cxmdk,handles.lists.Cxmdk,[],size(M,1));
%     clearInconsistentEntries(handles.lb_Cvmdk,handles.lists.Cvmdk,[],size(M,1));   
    
function lb_Dmdk_Callback(hObject, eventdata, handles)
% falls D=alpha*M+beta*K ausgew�hlt wurde, Felder f�r alpha und beta
% sichtbar machen, sonst unsichtbar, in UserData der Felder wird dann ggf
% wieder eine 1 eingetragen, wenn ein ung�ltiger Wert f�r alpha oder beta
% eingegeben wird
if get(hObject,'Value')==2
    set(handles.ed_a,'Visible','on')
    set(handles.ed_b,'Visible','on')
    set(handles.st_a,'Visible','on')
    set(handles.st_b,'Visible','on')
else
    set(handles.ed_a,'Visible','off')
    set(handles.ed_b,'Visible','off')
    set(handles.ed_a,'String','')
    set(handles.ed_b,'String','')
    set(handles.ed_a,'UserData',0)
    set(handles.ed_b,'UserData',0)
    set(handles.st_a,'Visible','off')
    set(handles.st_b,'Visible','off')
end


%--------------------------------------------------------------------------
%                            BOTTOM MENUE
%--------------------------------------------------------------------------

function ed_sysname_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
    % check variable names
    try
        isvalidvarname(hObject);
    end

function pb_create_Callback(hObject, eventdata, handles)
set(handles.figure1,'Pointer','watch');
k=[];
% �berpr�fen, ob �berall eine Matrix ausgew�hlt wurde (muss, wenn in der lb
% Eintr�ge vorhanden sind)
if get(handles.rb_ABCD,'Value')==1 && isempty(get(handles.lb_A,'String'))
    set(handles.figure1,'Pointer','arrow');
    errordlg('No A selected','Error Dialog','modal')
    uiwait
    return
end
if get(handles.rb_MDK,'Value')==1 && isempty(get(handles.lb_Mmdk,'String'))
    set(handles.figure1,'Pointer','arrow');
    errordlg('No M selected','Error Dialog','modal')
    uiwait
    return
end
if ~isvarname(get(handles.ed_sysname,'String'))
    set(handles.figure1,'Pointer','arrow');
    errordlg('Please correct System Name first','Error Dialog','modal')
    uiwait
    return
end
sysname=get(handles.ed_sysname,'String');
if exist_in_base_ws(sysname)~=0
%wenn es schon Variable im workspce mit diesem Namen gibt
          s=sprintf('%s already exists in base workspace. Do you want to overwrite it?',sysname);
          k=stqd('String',s,'Title','Question Dialog');
end
if ~isempty(k) && strcmp(k,'No')
    set(handles.figure1,'Pointer','arrow');
    return
end
% set(handles.figure1,'Pointer','watch');
if any(get(handles.ed_a,'UserData')==1) || any(get(handles.ed_b,'UserData')==1)
    set(handles.figure1,'Pointer','arrow');
    errordlg('Please correct alpha and beta first','Error Dialog','modal')
    uiwait
    return
end

if get(handles.rb_ABCD,'Value')==1
    
    % State-Space-Modell
    
    try
        A = get_matrix_from_ws(handles.lb_A);
        B = get_matrix_from_ws(handles.lb_B);
        C = get_matrix_from_ws(handles.lb_C);
        D = get_matrix_from_ws(handles.lb_D);
        E = get_matrix_from_ws(handles.lb_E);
        assignin('base',sysname,sss(A,B,C,D,E));
    catch ex
        errordlg(ex.message,'Error Dialog','modal')
        set(handles.figure1,'Pointer','arrow');
        return
    end
    set(handles.figure1,'Pointer','arrow');

    k=stqd('String','System created! Do you want to create another system?','Title','Question Dialog');
    if ~isempty(k) && strcmp(k,'No')
        figure1_DeleteFcn([], [], handles)
        return
    end
    
elseif get(handles.rb_MDK,'Value')==1
    
    % Second Order Systems
    
    try
        
        %Inform the User that the system will be transformed to first order
        
        dialog=stqd('String','The system will be converted from 2nd order form to 1st order form. Do you want to continue?','Title','Question Dialog');
        
        if ~isempty(dialog) && ~strcmp(dialog,'No')
        
            %Matrices M, K and B

            M = get_matrix_from_ws(handles.lb_Mmdk);
            K = get_matrix_from_ws(handles.lb_Kmdk);
            B = get_matrix_from_ws(handles.lb_Bmdk);

            %Matrix D

            Dsel = get(handles.lb_Dmdk, 'String');
            Dsel = Dsel{get(handles.lb_Dmdk, 'Value')};
            if strcmp(Dsel, 'alpha*M+beta*K')
                alpha=str2double(get(handles.ed_a,'String'));
                beta=str2double(get(handles.ed_b,'String'));
                if isnan(alpha) || isnan(beta)
                   error('Please correct the values for alpha and beta first!') 
                end
                D = alpha*M + beta*K;
            elseif strcmp(Dsel, 'Zero')
                D = [];
            else
                D = get_matrix_from_ws(handles.lb_Dmdk);
            end

            %Matrix F

            Fsel = get(handles.pu_Fmdk, 'String');
            Fsel = Fsel{get(handles.pu_Fmdk, 'Value')};
            
            if strcmp(Fsel,'I')
                Opts.transf2nd = 'I'; 
            elseif strcmp(Fsel,'K')
                Opts.transf2nd = 'K';
            else
                Opts.transf2nd = '-K';
            end

            %Matrix Cv

            Cvsel = get(handles.lb_Cvmdk, 'String');
            Cvsel = Cvsel{get(handles.lb_Cvmdk, 'Value')};
            if strcmp(Cvsel, 'Zero')
                Cv = [];
            else
                Cv = get_matrix_from_ws(handles.lb_Cvmdk);
            end

            %Matrix Cx

            Cxsel = get(handles.lb_Cxmdk, 'String');
            Cxsel = Cxsel{get(handles.lb_Cxmdk, 'Value')};
            if strcmp(Cxsel, 'Zero')
                Cx = [];
            else
                Cx = get_matrix_from_ws(handles.lb_Cxmdk);
            end

            %Conversion to first Order

            sys = second2first(M, D, K, B, Cx, Cv, Opts);
            assignin('base',sysname,sys);
            
            %User request
            
            set(handles.figure1,'Pointer','arrow');
            k=stqd('String','System created! Do you want to create another system?','Title','Question Dialog');% Dialogbox wir aufgerufen
            if isempty(k) || strcmp(k,'No') 
                figure1_DeleteFcn([], [], handles)
                return
            end
        
        end
        
    catch ex
        errordlg(ex.message,'Error Dialog','modal')
        set(handles.figure1,'Pointer','arrow');
        return
    end
    
    set(handles.figure1,'Pointer','arrow');

end
% MORLAB_GUI aktualisieren
try
    set(handles.inputGUI.lb_systems,'Value',[])
    set(handles.inputGUI.lb_systems,'String',systemsInWorkspace)
    set(handles.inputGUI.lb_matrixes,'Value',[])
    set(handles.inputGUI.lb_matrixes,'String',matricesInWorkspace)
catch exception
    if strcmp(exception.identifier,'MATLAB:class:InvalidHandle')
        % GUI wurde geschlossen, handles sind nicht mehr g�ltig
    else
        throw(exception)
    end
end



%--------------------------------------------------------------------------
%                          UTILITY FUNCTIONS
%--------------------------------------------------------------------------

function preselection(hObject,x)
% trifft Vorauswahl in der Create Systems GUI
% letzte �nderung 04.05.2010

%alle "zeilen" im String finden, in denen x voekommt
b=strfind(get(hObject,'String'),x);
%wenn es einen passenden eintrag in der datei gibt ist der genau x und
%kleiner als 6, weil nicht [ws] savor steht
for i=1:length(b)
    if length(b{i})==1 && b{i}<6
        set(hObject,'Value',i)
        return
    end
end
%wenn es keinen aus der datei gab, dann auch workspace eintr�ge
%ber�cksichtigen
for i=1:length(b)
    if length(b{i})==1
        set(hObject,'Value',i)
        return
    end
end
%wenn es auch im workspace nichts gibt, einen nehmen, in dem x vorkommt
%(z.B 'A_A' f�r A)
for i=1:length(b)
    if length(b{i})>1
        set(hObject,'Value',i)
        return
    end
end
%sonst auf den ersten eintrag setzten, dort stehen Einheitsmatrix oder Nullmatrix
set(hObject,'Value',1)

function i = exist_in_base_ws(varname)
% returns 1 if x exists in base workspace, otherwise 0
i=evalin('base',sprintf('exist(''%s'')',varname));

function mat = get_matrix_from_ws(handle)
% returns matrix selected in control handle

k=get(handle,'String');
l=regexp(k(get(handle,'Value')), '\ ', 'split'); % split string
if strcmp(l{1}{1,1},'[ws]')
    %variable im 'base' workspace
    varname=l{1}{1,2};
    mat=evalin('base',varname);
else
    %variable in load_var_dont_destroy
    varname=l{1}{1,1};
    try
        mat=evalin('base',sprintf('load_var_dont_destroy.%s',varname));
    catch
        mat=[];
    end
end

function list_matrices(handles)
% schreibt f�r die "create modela" GUI alle Matrizen aus dem Workspace/aus
% der Datei in die entsprechenden Listen der "create model" GUI

%letzte �nderung 27.12.2010

% x: Matrizen aus dem workspace
% x2: quadratische Matrizen aus dem workspace
% x3: Matrizen aus der geladenen Datei
% x4: quadratische Matrizen aus der geladenen Datei
if get(handles.cb_include,'Value')==1
    %variablen aus dem workspace dazu
    s=evalin('base', 'whos');
    % preallocate memory
    x=cell(length(s),1); %alle, auch quadtratische
    for i=1:length(s) 
        r=s(i).size;
        if strcmp(s(i).class,'double')&& length(r)==2 && any(r)
            % save name
            x{i}=sprintf('[ws] %s  [%ix%i]',s(i).name,s(i).size(1),s(i).size(2));
        end
    end
    % remove empty (non-system) entries
    x(cellfun(@isempty,x)) = [];
    % preallocate memory
    x2=cell(length(s),1); %quadratische, mindestens2x2
    for i=1:length(s)
        if strcmp(s(i).class,'double') && s(i).size(1)==s(i).size(2) && any(s(i).size>1) && length(s(i).size)==2
            % save name
            x2{i}=sprintf('[ws] %s  [%ix%i]',s(i).name,s(i).size(1),s(i).size(2));
        end
    end
    % remove empty (non-system) entries
    x2(cellfun(@isempty,x2)) = [];
else
    x={};
    x2={};
end
if get(handles.cb_file,'Value')==0
    %alte Version ~exist_in_base_ws('load_var_dont_destroy')
    %nur Variablen aus dem workspace
    x3={};
    x4={};
elseif exist_in_base_ws('load_var_dont_destroy')
% sonst auch Variablen aus Datei, die in load_var_dont_destroy gespeichert
% sind
    S=evalin('base','load_var_dont_destroy');
    f=fieldnames(S);
    x3=cell(length(f),1); %alle, auch quadtratische
    for i=1:length(f)
        current=sprintf('S.%s',f{i});
        r=size(eval(current));
        if strcmp(class(eval(current)),'double') && length(r)==2 && any(r)
            r=size(eval(current));
            % save name
            x3{i}=sprintf('%s  [%ix%i]',f{i},r(1),r(2));
        end
    end
    % remove empty (non-system) entries
    x3(cellfun(@isempty,x3)) = [];
    %s=whos;
    % preallocate memory
    x4=cell(length(f),1); %quadratische
    for i=1:length(f)
        current=sprintf('S.%s',f{i});
        r=size(eval(current));
        if strcmp(class(eval(current)),'double') && r(1)==r(2) && any(r>1) && length(r)==2
            % save name
            x4{i}=sprintf('%s  [%ix%i]',f{i},r(1),r(2));
        end
    end
    % remove empty (non-system) entries
    x4(cellfun(@isempty,x4)) = [];
% else
%     disp('Checkbox include file in list_matrixes.m �berlistet')
end
xquad={x2{:}, x4{:}};
xbel={x{:}, x3{:}};
if ~isempty(xquad)
    set(handles.pb_create,'Enable','on')
else
    set(handles.pb_create,'Enable','off')
end
%setzen entsprechend der Auswahl
set(handles.lb_A,'String', xquad) % alle quadratischen Matrizen
set(handles.lb_A,'Visible','on') 
preselection(handles.lb_A,'A') % sucht eine Matrix mit Namen A in der lb
set(handles.lb_B,'String', xbel)
set(handles.lb_B,'Visible','on')
preselection(handles.lb_B,'B')
set(handles.lb_C,'String', xbel)
set(handles.lb_C,'Visible','on')
preselection(handles.lb_C,'C')
set(handles.lb_D,'String', {'Zero',xbel{:}})
set(handles.lb_D,'Visible','on')
preselection(handles.lb_D,'D')
set(handles.lb_E,'String', {'identity matrix',xquad{:}})
set(handles.lb_E,'Visible','on')
preselection(handles.lb_E,'E')
set(handles.lb_Mmdk,'String',xquad)
preselection(handles.lb_Mmdk,'M')
set(handles.lb_Mmdk,'Visible','on')
set(handles.lb_Dmdk,'String', {'Zero','alpha*M+beta*K',xquad{:}})
preselection(handles.lb_Dmdk,'D')
set(handles.lb_Dmdk,'Visible','on')
set(handles.lb_Kmdk,'String', {'Zero',xquad{:}})
preselection(handles.lb_Kmdk,'K')
set(handles.lb_Kmdk,'Visible','on')
set(handles.lb_Bmdk,'String', xbel)
preselection(handles.lb_Bmdk,'B')
set(handles.lb_Bmdk,'Visible','on')
set(handles.lb_Cxmdk,'String', {'Zero',xbel{:}})
preselection(handles.lb_Cxmdk,'Cx')
set(handles.lb_Cxmdk,'Visible','on')
set(handles.lb_Cvmdk,'String', {'Zero',xbel{:}})
preselection(handles.lb_Cvmdk,'Cv')
set(handles.lb_Cvmdk,'Visible','on')
%set(handles.pu_Fmdk,'String', {'identity matrix',xquad{:}})
%preselection(handles.pu_Fmdk,'F')

%Select a matric named C if no Cx exists
if get(handles.lb_Cxmdk,'Value') == 1
   
    preselection(handles.lb_Cxmdk,'C');   
end

set(handles.pu_Fmdk,'Visible','on')
set(handles.pb_create,'Visible','on')
set(handles.ed_a,'Visible','off')
set(handles.ed_b,'Visible','off')
set(handles.st_a,'Visible','off')
set(handles.st_b,'Visible','off')

%Save the lists in the handles-structure

handles.lists.A = get(handles.lb_A,'String');
handles.lists.B = get(handles.lb_B,'String');
handles.lists.C = get(handles.lb_C,'String');
handles.lists.D = get(handles.lb_D,'String');
handles.lists.E = get(handles.lb_E,'String');
handles.lists.Mmdk = get(handles.lb_Mmdk,'String');
handles.lists.Dmdk = get(handles.lb_Dmdk,'String');
handles.lists.Kmdk = get(handles.lb_Kmdk,'String');
handles.lists.Bmdk = get(handles.lb_Bmdk,'String');
handles.lists.Cxmdk = get(handles.lb_Cxmdk,'String');
handles.lists.Cvmdk = get(handles.lb_Cvmdk,'String');


guidata(handles.lb_A,handles);

function [] = clearInconsistentEntries(hObject,list,size1,size2)
%Removes the entries which are not consistent with the given sizes from the
%list of a push-up-element 

    %Get the currently selected entry
    
    listOld = get(hObject,'String');
    sOld=listOld{get(hObject,'Value')};
    
    %Check for inconsistent entries
    
    for i = 1:length(list)
        if ~isempty(size1)      %Check dimension 1
            matSize = getMatrixSize(list{i,1},1);
            if matSize ~= size1 && matSize ~= 0
               list{i,1} = []; 
               continue;
            end
        end
        if ~isempty(size2)      %Check dimension 2
            matSize = getMatrixSize(list{i,1},2);
            if matSize ~= size2 && matSize ~= 0
               list{i,1} = []; 
            end
        end
    end
    
    %Remove the inconsistent entries from the list
    
    list = list(~cellfun(@isempty, list));
    
    %Select the same option as before if this entry is still there
    
    if ~isempty(list)

        indexOld=find(strcmp(sOld,list));

        if ~isempty(indexOld)
            set(hObject,'Value',indexOld);
        else
            set(hObject,'Value',1);
        end

    else
      list = {''};
      set(hObject,'Value',1);
    end
    
    %Write the modified list to the push-up-element
    
    set(hObject,'String',list)

function x = getMatrixSize(name,dim)
%Returns the size of the dimension dim of the matrix corresponding to the
%given name

    l=regexp(name, '\ ', 'split'); % split string
    if strcmp(l{1,1},'[ws]')
        %variable im 'base' workspace
        varname=l{1,2};
        mat=evalin('base',varname);
    else
        %variable in load_var_dont_destroy
        varname=l{1,1};
        try
            mat=evalin('base',sprintf('load_var_dont_destroy.%s',varname));
        catch
            mat=[];
        end
    end
    
    if dim == 1
       x = size(mat,1); 
    else
       x = size(mat,2);
    end

