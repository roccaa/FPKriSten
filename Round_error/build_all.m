
%% Build all the examples
cd Examples/

files_names = {'simple_test';'rigidBody1';'rigidBody2';'kepler0real';'kepler1';'kepler2';'sineTaylor';'sineOrder3';'sqroot';'himmilbeau'};

for i=1:size(files_names,1)
    disp(['Building ' files_names{i} '...']);
    eval([files_names{i} '_m']);
end

cd ..