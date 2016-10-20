tag = 'MAX'

files_names = {'simple_test';'rigidBody1';'rigidBody2';'kepler0real';'kepler1';'kepler2';'sineTaylor';'sineOrder3';'sqroot';'himmilbeau';'schwefel';'magnetism';'caprasse'};

result = cell(5 ,size(files_names,1)+1);
result{1,1} = 'Name  ';
result{2,1} = 'Roundofferror  ';
result{3,1} = 'Building time  ';
result{4,1} = 'Solving time  ';
result{5,1} = 'Total time  ';
for i=1:size(files_names,1)
    [F,I,J ,G ,n,d,k] = read_examples(files_names{i},tag);
    [ bound,build_time,solving_time ] = solve_examples( F,G,I,J,d,k,tag);
    result{1,i+1} = files_names{i};
    result{2,i+1} = bound*2^(-53);
    result{3,i+1} = build_time;
    result{4,i+1} = solving_time;
    result{5,i+1} = solving_time+build_time;
end
