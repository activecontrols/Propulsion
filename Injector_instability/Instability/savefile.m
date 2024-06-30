function savefile(dpf_pc,dpo_pc,dpf_pc1,dpf_pc2,dpo_pc1,dpo_pc2,f,h2)
% j starts on 1 so that the first row is zeros. This is because the header
% has to write on row 1 and there is no easy way to append the data or
% start the data on row 2. Row 1 is all zeros but the header writes over
% this. The same process is used for Excel format (but not necessary) to
% keep the code from rebuilding the Result matrix starting from a different
% row.
if get(h2,'Value')==1 % A value of 1 is Feed System On
    j=1;
    for k=1:length(dpf_pc1)
        if (dpf_pc1(k)>0 && isnan(dpf_pc1(k))==0 && isnan(dpo_pc1(k))==0)
            j=j+1;
            Result(j,1)=dpo_pc1(k);
            Result(j,2)=dpf_pc1(k);
            Result(j,3)=f(k);
        end
    end
    for k=length(dpf_pc2)+1:2*length(dpf_pc2)
        if (dpf_pc2(k-length(dpf_pc1))>0 && isnan(dpf_pc2(k-length(dpf_pc1)))==0 && isnan(dpo_pc2(k-length(dpf_pc1)))==0)
        j=j+1;
        Result(j,1)=dpo_pc2(k-length(dpf_pc1));
        Result(j,2)=dpf_pc2(k-length(dpf_pc1));
        Result(j,3)=f(k-length(dpf_pc1));
        end
    end
    if j<=65535 % then write in excel
        header={'DPo/Pc' 'DPf/Pc' 'freq (Hz)'};
        xlswrite('FeedSystemOn.xls',Result,'Sheet1','A1');
        xlswrite('FeedSystemOn.xls',header,'Sheet1','A1');
    else % otherwise write into a text file
        % Save logic for script below: NO easy way to save both string and
        % number data to a file. So save number data first starting on row 2.
        % Then open file without discarding contents and write the text
        % header on the first row.
        % Save data in column ascii format. Ignore NaN and 0 values.
        save FeedSystemOn.out Result -ascii
        % This opens file for reading and writing without discarding contents and
        % then saves a header to the file.
        fid=fopen('FeedSystemOn.out','r+');
        fprintf(fid,'%s %s %s',[' DPo/Pc ',' DPf/Pc ',' Freq(Hz) ']);
        fclose(fid);
    end
else % A value of 0 is Feed System Off
    % Save data in column ascii format. Ignore NaN and 0 values.
    j=1;
    for k=1:length(dpf_pc)
    if (dpf_pc(k)>0 && isnan(dpf_pc(k))==0 && isnan(dpo_pc(k))==0)
    j=j+1;
    Result(j,1)=dpo_pc(k);
    Result(j,2)=dpf_pc(k);
    Result(j,3)=f(k);
    end
    end
    if j<=65535 % then write in excel
    header={'DPo/Pc' 'DPf/Pc' 'freq (Hz)'};
    xlswrite('FeedSystemOff.xls',Result,'Sheet1','A1'); % first row is all zeros
    xlswrite('FeedSystemOff.xls',header,'Sheet1','A1');
    else % otherwise write into a text file
    save FeedSystemOff.out Result -ascii % save number data
    % This opens file for reading and writing without discarding contents and
    % then saves a header to the file.
    fid=fopen('FeedSystemOff.out','r+');
    fprintf(fid,'%s %s %s',[' DPo/Pc ',' DPf/Pc ',' Freq(Hz) ']);
    fclose(fid);
    end
end