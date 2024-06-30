% This function performs actions when the radiobuttons are selected
function selcbk(source,eventdata,s1)
global fileinp a1 a2 a3 a4
% Each radio button is tagged. This greys out the edit box that is not
% used for that particular radio button.
% disp(get(get(source,'SelectedObject'),'String'));
if str2num(get(eventdata.OldValue,'Tag'))==1
if fileinp==0
set(s1.edit22, 'BackgroundColor', 'white', 'Enable', 'on','string', a1);
set(s1.edit23, 'BackgroundColor', 'white', 'Enable', 'on','string', a2);
set(s1.edit24, 'BackgroundColor', 'white', 'Enable', 'on','string', a3);
set(s1.edit25, 'BackgroundColor', 'white', 'Enable', 'on','string', a4);
else
set(s1.edit22, 'BackgroundColor', 'white', 'Enable', 'off','string', a1);
set(s1.edit23, 'BackgroundColor', 'white', 'Enable', 'off','string', a2);
set(s1.edit24, 'BackgroundColor', 'white', 'Enable', 'off','string', a3);
set(s1.edit25, 'BackgroundColor', 'white', 'Enable', 'off','string', a4);
end
else
a1=get(s1.edit22, 'string');
a2=get(s1.edit23, 'string');
a3=get(s1.edit24, 'string');
a4=get(s1.edit25, 'string');
set(s1.edit22, 'BackgroundColor', 'white', 'Enable', 'off','string', '0');
set(s1.edit23, 'BackgroundColor', 'white', 'Enable', 'off','string', 'inf.');
set(s1.edit24, 'BackgroundColor', 'white', 'Enable', 'off','string', '0');
set(s1.edit25, 'BackgroundColor', 'white', 'Enable', 'off','string', 'inf.');
end
