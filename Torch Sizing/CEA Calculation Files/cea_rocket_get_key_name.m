function [key_name_new] = cea_rocket_get_key_name(key_name_old, is_lower)
%CEA_ROCKET_GET_KEY_NAME Summary of this function goes here
%   Detailed explanation goes here

if ~exist('is_lower', 'var') || is_lower
    key_name_old = lower(strtrim(key_name_old));
else
    key_name_old = strtrim(key_name_old);
end

key_space = strfind(key_name_old, ' ');
if isempty(key_space)
    key_space(1) = length(key_name_old)+1;
end

key_comma = strfind(key_name_old, ',');
if isempty(key_comma)
    key_comma(1) = length(key_name_old)+1;
end

key_name_new = key_name_old(1:min([key_space(1)-1, key_comma(1)-1, ...
    length(key_name_old)]));

end