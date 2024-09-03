function [ErrEmu] = FindErEMU(memu)
%load('debug');
emuSiz = numel(memu{1,1}.emu);
netWsiz = numel(memu{1,1}.A);
p = 0;
ErrEmu={};
for i = netWsiz:-1:1
    name_{i,1} = str2var(strcat('x',num2str(i)));
    name_{i,2} = sum(name_{i,1},2);
    Rowss = size(name_{i,1},1);
    for k = 1:Rowss
        if name_{i,2}(k,1) > 1.001 || name_{i,2}(k,1) < 0.999
            for j = emuSiz:-1:1
                if memu{1,1}.emu{j,1}.network == i && memu{1,1}.emu{j,1}.index == k
                    p = p+1;
                    ErrEmu{p,1} = memu{1,1}.emu{j,1}.id;
                    ErrEmu{p,2} = memu{1,1}.emu{j,1}.cfg;
                    ErrEmu{p,3} = i;
                    ErrEmu{p,4} = k;
                    ErrEmu{p,5} = name_{i,2}(k,1);
                end
            end
        end
    end
end
end

