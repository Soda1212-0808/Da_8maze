
newdata.date_2024_05_30=[animal_match_table.date_2024_05_30 table2cell(animal_path.date_2024_05_30(cell2mat(animal_match_table.date_2024_05_30(:,3))+1,:)  )];


% 获取animal_match_table的所有字段名（即所有日期）
dateFields = fieldnames(animal_match_table);

% 定义一个匿名函数来对每个日期字段进行操作
operateOnDate = @(dateField) ...
    [animal_match_table.(dateField) ...
    table2cell(animal_path.(dateField)(cell2mat(animal_match_table.(dateField)(:,3)) + 1, :))];

% 使用cellfun对所有日期字段进行操作，并将结果存入newdata结构体中
newdata = cell2struct(cellfun(operateOnDate, dateFields, 'UniformOutput', false), dateFields);
