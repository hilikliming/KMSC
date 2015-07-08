function [ PondY ] = getPondTemplates( pondDir,targs )
cd(pondDir);

x = what;
x = x.mat;
for o = 1:length(targs)
    for i = 1:length(x)
        if(strfind(char(x(i),targs(o))
            ob = open(x(i));
            PondY(o).D = ob.acp;
        end
    end
end

% Resampling Templates to 310 frequency bins for 0-31 kHz....



end

