clear;
close all;

inputfilepath='.\datasets\';
file_row='GeneCoexpWang_genelist.txt';
file_col='GeneCoexpWang_genelist.txt';
fid_row=fopen(strcat(inputfilepath, file_row),'r');
fid_col=fopen(strcat(inputfilepath, file_col));
rownamelist=textscan(fid_row, '%s');
colnamelist=textscan(fid_col, '%s', 'delimiter', '\n');
fclose(fid_row);
fclose(fid_col);

patternfilepath='.\survival\';
file_patternid='Wang_0601_07_2nd_S.survival';
fid_pattern=fopen(strcat(patternfilepath, file_patternid));
line=fgetl(fid_pattern);
patterns_rowinfo=[];
patterns_colinfo=[];
patterns_pvalueinfo=[];
while ischar(line);
    [sid eid]=regexpi(line, ';');% Assume there are two semicolons each line
    line_left=line(1:sid(1)-1);
    line_right=line(sid(1)+1:sid(2)-1);
    line_pvalue=line(sid(2)+1:end);
    row_genes=textscan(line_left, '%s');
    col_genes=textscan(line_right,'%s');
    patterns_rowinfo=[patterns_rowinfo, row_genes];
    patterns_colinfo=[patterns_colinfo, col_genes];
    patterns_pvalueinfo=[patterns_pvalueinfo, {line_pvalue}];
    line=fgetl(fid_pattern);
end
fclose(fid_pattern);

file_encoded='Wang_0601_07_2nd_S.encoded';
fid_encoded=fopen(strcat(patternfilepath, file_encoded), 'wt');
for k=1:length(patterns_rowinfo)
    for i=1:length(patterns_rowinfo{k})
        inds=strmatch(patterns_rowinfo{k}(i), rownamelist{1}, 'exact');
        if(length(inds)>1)
            fprintf('one gene correponds to more than one row.\n');
        end
        fprintf(fid_encoded, '%d ', inds(1));
    end
    fprintf(fid_encoded, ';');
    for i=1:length(patterns_colinfo{k})
        inds=strmatch(patterns_colinfo{k}(i), colnamelist{1}, 'exact');
        if(length(inds)>1)
            fprintf('one gene correponds to more than one col.\n');
        end
        fprintf(fid_encoded, '%d ', inds(1));
    end
    fprintf(fid_encoded, ';');
    fprintf(fid_encoded, '%s\n', patterns_pvalueinfo{k});
end
fclose(fid_encoded);