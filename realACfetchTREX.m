function [ Y, t_Y, r_Y ] = realACfetchTREX( realTarg,stops,eps,rngs )
home     = cd;
targetID = 1;
Y   = []; 
t_Y = [];
r_Y = [];
upperF  = 31e3; % Chosen based on SERDP MR-1665-FR Final Report
f_s     = 100e3;% determined experimentally...
sig_N   = 310;  % (0-31 kHz)
%eps      = 55; %Experimentally determined threholding value for grabbing important aspects
%Dclutter = [];
for tag = realTarg
    for rng = rngs(targetID,:)
        if(rng>0)
            cd(['F:\TREX13\TARGET_DATA\',char(tag),'\',num2str(rng),'m\sand\proud',]);
            here = cd; % Where the object run .mat files are located
            x = what; x = x.mat;
            burysign = 1;
            if(length(x)<1)
                cd(['F:\TREX13\TARGET_DATA\',char(tag),'\',num2str(rng),'m\sand\part_buried',]);
                here = cd;
                x = what; x = x.mat;
                burysign = -1;
            end
            % Various rotations of object are captured in each run
            for run = 1:length(x)
                if(~isempty(strfind(char(x(run)),'_norm'))) % if it's well conditioned type
                    ob=open(char(x(run))); % char() just converts the cell to string
                    cd(home);
                    %load('APL_VLA_1_30kHz.mat');
                    %xmit = decimate(data_V,10);
                    pings = ob.new_data';%PulseCompress(ob.new_data',xmit,f_s,0);%
                    AC = extractAC(pings,eps,sig_N,upperF,f_s,stops);
                    Y   = [Y, AC];
                    t_Y = [t_Y', burysign*targetID*ones(size(AC,2),1)']';
                    r_Y = [r_Y', rng*ones(size(AC,2),1)']';
                    cd(here);
                end
            end
        end
    end
    targetID = targetID + 1;  
end
cd(home);
end
