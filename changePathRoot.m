function [ NewdirMapDB ] = changePathRoot( new_root,dirMapDB )
%This is a simple script to change the directory map (used by
%getTrainingSamples) to have the root path changed to something other than
%the original Map directs to (e.g. change E:\DBFRM\ENV_1\aluxo\10_m\0_deg
%... to C:\Users\JonGuy\Matlab\DBFRM\ENV_1\aluxo\10_m\0_deg
NewdirMapDB = dirMapDB;
for e = 1:size(dirMapDB,1)
    for o = 1:size(dirMapDB,2)
        for r = 1:size(dirMapDB,3)
            for rot = 1:size(dirMapDB,4)
                dir = char(dirMapDB(e,o,r,rot));
                root = strfind(char(dir),'DBFRM\');
                dir = dir(root:end);
                NewdirMapDB(e,o,r,rot) = cellstr([new_root,'\',dir]);
            end
        end
    end
end


end

