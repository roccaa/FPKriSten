
%% Build all the examples
cd Examples/

files_names = {'simple_test';'rigidBody1';'rigidBody2';'kepler0real';'kepler1';...
               'kepler2';'sineTaylor';'sineOrder3';'sqroot';'himmilbeau';'schwefel';...
               'magnetism';'caprasse';'exemple_2_2_5';'exemple_2_2_10';'exemple_2_2_15';...
               'exemple_2_2_20';'exemple_2_5_2';'exemple_2_10_2';'exemple_5_2_2';'exemple_10_2_2';...
               'floudas2_6';'floudas3_3';'floudas3_4';'floudas4_6';'floudas4_7'};
%'floudas2_6';
for i=1:size(files_names,1)
    disp(['Building ' files_names{i} '...']);
    eval([files_names{i} '_m']);
end

cd ..