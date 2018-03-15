tag = 'MAX'

files_names = {'simple_test';'rigidBody1';'rigidBody2';'kepler0real';'kepler1';...
               'kepler2';'sineTaylor';'sineOrder3';'sqroot';'himmilbeau';'schwefel';...
               'magnetism';'caprasse';'exemple_2_2_5';'exemple_2_2_10';'exemple_2_2_15';...
               'exemple_2_2_20';'exemple_2_5_2';'exemple_2_10_2';'exemple_5_2_2';'exemple_10_2_2';...
               'floudas3_3';'floudas3_4';'floudas4_6';'floudas4_7'};
%'floudas2_6';
result = cell(5 ,size(files_names,1)+1);
result{1,1} = 'Name  ';
result{2,1} = 'Roundofferror  ';
result{3,1} = 'Building time  ';
result{4,1} = 'Solving time  ';
result{5,1} = 'Total time  ';
nb_times = 5;
for i=1:size(files_names,1)
    [F,I,J ,G ,n,d,k] = read_examples(files_names{i},tag);
    tb = 0;
    ts = 0;
    tt = 0;
    for j=1:nb_times
        [bound,build_time,solving_time ] = solve_examples( F,G,I,J,d,k,tag);
        tb = tb + build_time;
        ts = ts + solving_time;
        tt = tt + solving_time+build_time;
    end
        result{1,i+1} = files_names{i};
        result{2,i+1} = bound*2^(-53);
        result{3,i+1} = tb/5;
        result{4,i+1} = ts/5;
        result{5,i+1} = tt/5;
end
