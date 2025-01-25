clc; 
clear all; 
close all;
[num, text, raw] = xlsread('ipadata.xlsx'); 
size = size(text); 
Datamatrix = zeros(size(1),14);
for i = 1:size(1)
    string = text(i,1); 
    splitString = split(string);
    numbersarray = transpose(str2double(splitString));
    Datamatrix(i,1:length(numbersarray)) = numbersarray;
end
%Cells = num2cell(num)
%writecell(Cells,'ipadata.xlsx','Sheet',3,'Range','A2:G38')
a = "done"