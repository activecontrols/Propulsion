function [new_file_name] = rep_drive_letter(file_name)

% Get the drive mappings.
[~, cmd_out] = system('net use');

% Iterate through each line of the command output.
lines = splitlines(cmd_out);
for i = 1:length(lines)
    % The line is a mapping if it starts with "OK".
    if startsWith(lines{i}, 'OK')
        i_colon = strfind(lines{i}, ':');
        letter = lines{i}((i_colon-1):i_colon);
        UNC = strtrim(lines{i}((i_colon+1):end));
        if exist('file_name', 'var')
            if startsWith(file_name, UNC)
                new_file_name = [letter, file_name((length(UNC)+1):end)];
                return;
            end
        end
    end
end

new_file_name = file_name;

end