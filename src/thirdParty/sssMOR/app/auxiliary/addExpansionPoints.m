function x=addExpensionPoints(x,eventdata)
% Diese Funktion soll die multiple expansion points Liste bearbeiten
% 1. Fall: eventdata==[]
%     �berpr�fen ob f�r alle punkte auch momente angegeben sind, wenn nicht
%         wird eine 1 eingetragen
%     gleiche punkte zusammenfassen und momente-zahl entsprechen erh�hen
% 2. Fall: eventdata enth�lt folgende Eintr�ge:
%     %eventdata  structure with the following fields (see UITABLE)
%     %	Indices: row and column indices of the cell(s) edited
%     %	PreviousData: previous data for the cell(s) edited
%     %	EditData: string(s) entered by the user
%     %	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%     %	Error: error string when failed to convert EditData to appropriate
%     %	value for Data
%     Funktion wurde vom edit callback der Liste aufgerufen
%     1. wenn ein Eintrag gel�scht werden soll, l�schen und ggf den konjugiert
%         komplexen Partner mitl�schen
%     2. Punkte die ohne eventdata auch erledigt werden
%     3. Konjugiert komplexe Paare erg�nzen, bzw Momente der beiden anpassen

% letzte �nderung 11.08.2010

if isempty(eventdata)%gleiche zusammenfassen und matching moments auff�llen, 
    %konj komplexe Paare erg�nzen/momente anpassen
    for i=1:size(x,1)
        if ~isempty(x{i,1}) && isempty(x{i,2}) %wenn nur expt ausgef�llt wurde 
            %aber noch kein matching mom, auf eins setzen
            x{i,2}=1;
        end
    end
    for i=1:size(x,1) %gleiche zusammenfassen
        for j=i+1:size(x,1)
            if x{i,1}==x{j,1}
                x{i,2}=x{i,2}+x{j,2};
                x{j,1}=[];x{j,2}=[];
            end
        end
    end
    for i=1:size(x,1)
        k=[];
        if imag(x{i,1})~=0 %konj kompl paar einf�gen
            for j=1:size(x,1)
                if ~isempty(x{j,1}) && conj(x{i,1})==x{j,1}
                    l=0;
                    x{i,2}=x{j,2};
                    break
                end
                l=1;
            end
            if l==1 %konj kompl nicht vorhanden            
                for j=1:size(x,1) %n�chste freie Zeile suchen
                    if isempty(x{j,1})
                        k=j;
                        break
                    end
                end
                if isempty(k) %wenn keine frei, neue erforlderlich
                    k=size(x,1)+1;
                end
                x{k,1}=conj(x{i,1});
                x{k,2}=x{i,2};
            end
        end
    end
    return
end


if isempty(eventdata.EditData)%Zeile soll gel�scht werden
    if eventdata.Indices(2)==1
        x_alt=eventdata.PreviousData;
    else
        x_alt=x{eventdata.Indices(1),1};
    end
    x{eventdata.Indices(1),1}=[];
    x{eventdata.Indices(1),2}=[];
    if imag(x_alt)~=0%konj kompl suchen und l�schen
        for i=1:size(x,1)
            if x{i,1}==conj(x_alt)
                x{i,1}=[];
                x{i,2}=[];
                break
            end
        end
    end
end
for i=1:size(x,1)
    if ~isempty(x{i,1}) && isempty(x{i,2}) %wenn nur expt ausgef�llt wurde 
        %aber noch kein matching mom, auf eins setzen
        x{i,2}=1;
    end
end
for i=1:size(x,1) %gleiche zusammenfassen
    for j=i+1:size(x,1)
        if x{i,1}==x{j,1}
            x{i,2}=x{i,2}+x{j,2};
            x{j,1}=[];x{j,2}=[];
        end
    end
end
for i=1:size(x,1)
    k=[];
    if imag(x{i,1})~=0 %konj kompl paar einf�gen
        for j=1:size(x,1)
            if ~isempty(x{j,1}) && conj(x{i,1})==x{j,1}
                l=0;%konj kompl schon vorhanden
                if eventdata.Indices(2)==2 %match mom wurden ver�ndert
                    if eventdata.Indices(1)==i
                        x{j,2}=x{i,2};
                    elseif eventdata.Indices(1)==j
                        x{i,2}=x{j,2};
                    else
                        x{j,2}=x{i,2};
                    end
                else %ein schon vorhandenes wurde erneut dazu geschrieben
                    if x{j,2}>x{i,2}
                        x{i,2}=x{j,2};
                    else
                        x{j,2}=x{i,2};
                    end
                end
                break
            end
            l=1;
        end
        if l==1 %konj kompl nicht vorhanden            
            for j=1:size(x,1) %n�chste freie Zeile suchen
                if isempty(x{j,1})
                    k=j;
                    break
                end
            end
            if isempty(k) %wenn keine frei, neue erforlderlich
                k=size(x,1)+1;
            end
            x{k,1}=conj(x{i,1});
            x{k,2}=x{i,2};
        end
    end
end
