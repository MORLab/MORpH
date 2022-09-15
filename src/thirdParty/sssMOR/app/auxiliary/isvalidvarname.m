function isvalidvarname(hObject,varargin)
% �berpr�ft den String des �bergebenen Objekts auf G�ltigkeit des
% Variablennamens und ob schon eine Variable mit diesem Namen im workspace
% liegt
% �ber varargin{1}=='NoWarning' kann warnung unterdr�ckt werden

% letzte �nderung 22.10.2010
if nargin<2
    varargin{1}=[];
end
if ~isvarname(get(hObject,'String'))
    errordlg('Has to be a valid variable Name','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
    return
elseif existInBaseWs(get(hObject,'String'))==1  && ~strcmp(varargin{1},'NoWarning')
    warndlg('Variable with this name already exists in workspace','Warning','modal')
    uiwait
end
set(hObject,'UserData',0)