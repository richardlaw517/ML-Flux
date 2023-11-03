function [eq,textrep] = xlsreadeq(xlsname,model,net,xch)
%% read excel file containing equality constraints
% Dependencies: consfree2full.m
% input: xlsname (excel file name) containing equalities (==); one side
% must contain only reaction names and the other side, only numbers; model
% (MATLAB structure containing metabolic model); net (vector of free net
% fluxes); xch (vector of free exchange fluxes)
% output: eq (structure containing A and b (such that A*free_fluxes==b));
% textrep (text representation; currently not used)

[~,~,N,all_free,~,rm_cons]=consfree2full(model.kernel_net,model.kernel_xch,model.fluxes,net,xch);
rxns=model.rxns;
rxns_l=length(rxns);
colN=size(N,2);

if ~isempty(xlsname)
    sheets={'net_eq' 'xch_eq'};
    for s=1:2
        if ispc
            [~,~,alldata]=xlsread(xlsname,sheets{s});
        elseif isunix
            warning('off','MATLAB:xlsread:Mode')
            [~,~,alldata]=xlsread(xlsname,sheets{s},'','basic');
        end
        c1=s-1;
        alldata(~(cellfun(@ischar,alldata(:,1))|cellfun(@ischar,alldata(:,3))),:)=[];
        alldata_l=size(alldata,1);
        A2.(sheets{s})=zeros(alldata_l,size(N,2));
        b2.(sheets{s})=zeros(alldata_l,1);
        textrep2.(sheets{s})=cell(alldata_l,1);
        for i=1:alldata_l
            if ischar(alldata{i,1}) && ~ischar(alldata{i,3})
%                 % remove spaces and punctuation
%                 alldata{i,1}(ismember(alldata{i,1},' ,:;!')) = [];
%                 textrep2.(sheets{s}){i}=alldata{i,1};
%                 psub=find(alldata{i,1}=='-');
%                 nsub=length(psub);
%                 if nsub==0
%                     A2.(sheets{s})(i,:)= N([false(1,c1*rxns_l) strcmp(alldata{i,1},rxns)],:);
%                     b2.(sheets{s})(i)= alldata{i,3};
%                 else % when there is 1 or more - operators
%                     A2.(sheets{s})(i,:)= N([false(1,c1*rxns_l) strcmp(alldata{i,1}(1:psub(1)-1),rxns)],:);
%                     b2.(sheets{s})(i)= alldata{i,3};
%                     for j=1:nsub-1
%                         A2.(sheets{s})(i,:)=A2.(sheets{s})(i,:)-N([false(1,c1*rxns_l) strcmp(alldata{i,1}(psub(j)+1:psub(j+1)-1),rxns)],:);
%                     end
%                     A2.(sheets{s})(i,:)=A2.(sheets{s})(i,:)-N([false(1,c1*rxns_l) strcmp(alldata{i,1}(psub(nsub)+1:end),rxns)],:);
%                 end
                % remove spaces and punctuation
                alldata{i,1}(ismember(alldata{i,1},' ,:;!')) = [];
                textrep2.(sheets{s}){i}=[alldata{i,1} '==' num2str(alldata{i,3})];
                if alldata{i,1}(1)~='+' && alldata{i,1}(1)~='-'
                    alldata{i,1}=['+' alldata{i,1}];
                end
                psub1=find(alldata{i,1}=='-' | alldata{i,1}=='+' | alldata{i,1}=='/' | alldata{i,1}=='*');
                nsub1=length(psub1);
                psub1=[psub1 length(alldata{i,1})+1];
                
                b2.(sheets{s})(i)=alldata{i,3};
                L=zeros(1,colN);
                
                for j=1:nsub1
                    switch alldata{i,1}(psub1(j))
                        case '+'
                            if j<nsub1 && (alldata{i,1}(psub1(j+1))=='*' || alldata{i,1}(psub1(j+1))=='/')
                                if sum(strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns))
                                    temp=N([false(1,c1*rxns_l) strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1)))
                                    temp=str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1));
                                end
                            else
                                if sum(strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns))
                                    L=L+N([false(1,c1*rxns_l) strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1)))
                                    b2.(sheets{s})(i)=b2.(sheets{s})(i)-str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1));
                                end
                            end
                        case '-'
                            if j<nsub1 && (alldata{i,1}(psub1(j+1))=='*' || alldata{i,1}(psub1(j+1))=='/')
                                if sum(strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns))
                                    temp=-N([false(1,c1*rxns_l) strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1)))
                                    temp=-str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1));
                                end
                            else
                                if sum(strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns))
                                    L=L-N([false(1,c1*rxns_l) strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1)))
                                    b2.(sheets{s})(i)=b2.(sheets{s})(i)+str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1));
                                end
                            end
                        case '*'
                            if sum(strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns)) && isscalar(temp)
                                L=L+temp*N([false(1,c1*rxns_l) strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns)],:);
                            elseif ~isnan(str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1)))
                                L=L+temp*str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1));
                            else
                                warning([sheets{s} ' ' num2str(i) ': nonlinear constraint'])
                            end
                        case '/'
                            if sum(strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns))
                                warning([sheets{s} ' ' num2str(i) ': nonlinear constraint'])
                            elseif ~isnan(str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1)))
                                L=L+temp/str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1));
                            end
                        otherwise
                            warning([sheets{s} ' ' num2str(i) ': unexpected operation'])
                    end
                end
                A2.(sheets{s})(i,:)=L;
                
            elseif ~ischar(alldata{i,1}) && ischar(alldata{i,3})
%                 % remove spaces and punctuation
%                 alldata{i,3}(ismember(alldata{i,3},' ,:;!')) = [];
%                 textrep2.(sheets{s}){i}=alldata{i,3};
%                 psub=find(alldata{i,3}=='-');
%                 nsub=length(psub);
%                 if nsub==0
%                     A2.(sheets{s})(i,:)=-N([false(1,c1*rxns_l) strcmp(alldata{i,3},rxns)],:);
%                     b2.(sheets{s})(i)=-alldata{i,1};
%                 else % when there is 1 or more - operators
%                     A2.(sheets{s})(i,:)=-N([false(1,c1*rxns_l) strcmp(alldata{i,3}(1:psub(1)-1),rxns)],:);
%                     b2.(sheets{s})(i)=-alldata{i,1};
%                     for j=1:nsub-1
%                         A2.(sheets{s})(i,:)=A2.(sheets{s})(i,:)+N([false(1,c1*rxns_l) strcmp(alldata{i,3}(psub(j)+1:psub(j+1)-1),rxns)],:);
%                     end
%                     A2.(sheets{s})(i,:)=A2.(sheets{s})(i,:)+N([false(1,c1*rxns_l) strcmp(alldata{i,3}(psub(nsub)+1:end),rxns)],:);
%                 end
                alldata{i,3}(ismember(alldata{i,3},' ,:;!')) = [];
                textrep2.(sheets{s}){i}=[num2str(alldata{i,1}) '==' alldata{i,3}];
                if alldata{i,3}(1)~='+' && alldata{i,3}(1)~='-'
                    alldata{i,3}=['+' alldata{i,3}];
                end
                psub2=find(alldata{i,3}=='-' | alldata{i,3}=='+' | alldata{i,3}=='/' | alldata{i,3}=='*');
                nsub2=length(psub2);
                psub2=[psub2 length(alldata{i,3})+1];
                
                b2.(sheets{s})(i)=-alldata{i,1};
                R=zeros(1,colN);
                
                for j=1:nsub2
                    switch alldata{i,3}(psub2(j))
                        case '+'
                            if j<nsub2 && (alldata{i,3}(psub2(j+1))=='*' || alldata{i,3}(psub2(j+1))=='/')
                                if sum(strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns))
                                    temp=N([false(1,c1*rxns_l) strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1)))
                                    temp=str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1));
                                end
                            else
                                if sum(strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns))
                                    R=R+N([false(1,c1*rxns_l) strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1)))
                                    b2.(sheets{s})(i)=b2.(sheets{s})(i)-str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1));
                                end
                            end
                        case '-'
                            if j<nsub2 && (alldata{i,3}(psub2(j+1))=='*' || alldata{i,3}(psub2(j+1))=='/')
                                if sum(strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns))
                                    temp=-N([false(1,c1*rxns_l) strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1)))
                                    temp=-str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1));
                                end
                            else
                                if sum(strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns))
                                    R=R-N([false(1,c1*rxns_l) strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1)))
                                    b2.(sheets{s})(i)=b2.(sheets{s})(i)+str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1));
                                end
                            end
                        case '*'
                            if sum(strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns)) && isscalar(temp)
                                R=R+temp*N([false(1,c1*rxns_l) strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns)],:);
                            elseif ~isnan(str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1)))
                                R=R+temp*str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1));
                            else
                                warning([sheets{s} ' ' num2str(i) ': nonlinear constraint'])
                            end
                        case '/'
                            if sum(strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns))
                                warning([sheets{s} ' ' num2str(i) ': nonlinear constraint'])
                            elseif ~isnan(str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1)))
                                R=R+temp/str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1));
                            end
                        otherwise
                            warning([sheets{s} ' ' num2str(i) ': unexpected operation'])
                    end
                end
                A2.(sheets{s})(i,:)=-R;
                
            elseif ischar(alldata{i,1}) && ischar(alldata{i,3})
                % remove spaces and punctuation
                alldata{i,1}(ismember(alldata{i,1},' ,:;!')) = [];
                alldata{i,3}(ismember(alldata{i,3},' ,:;!')) = [];
                textrep2.(sheets{s}){i}=[alldata{i,1} '==' alldata{i,3}];
                
                if alldata{i,1}(1)~='+' && alldata{i,1}(1)~='-'
                    alldata{i,1}=['+' alldata{i,1}];
                end
                psub1=find(alldata{i,1}=='-' | alldata{i,1}=='+' | alldata{i,1}=='/' | alldata{i,1}=='*');
                nsub1=length(psub1);
                psub1=[psub1 length(alldata{i,1})+1];
                
                if alldata{i,3}(1)~='+' && alldata{i,3}(1)~='-'
                    alldata{i,3}=['+' alldata{i,3}];
                end
                psub2=find(alldata{i,3}=='-' | alldata{i,3}=='+' | alldata{i,3}=='/' | alldata{i,3}=='*');
                nsub2=length(psub2);
                psub2=[psub2 length(alldata{i,3})+1];
                
                b2.(sheets{s})(i)=0;
                L=zeros(1,colN);
                R=L;

                for j=1:nsub1
                    switch alldata{i,1}(psub1(j))
                        case '+'
                            if j<nsub1 && (alldata{i,1}(psub1(j+1))=='*' || alldata{i,1}(psub1(j+1))=='/')
                                if sum(strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns))
                                    temp=N([false(1,c1*rxns_l) strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1)))
                                    temp=str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1));
                                end
                            else
                                if sum(strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns))
                                    L=L+N([false(1,c1*rxns_l) strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1)))
                                    b2.(sheets{s})(i)=b2.(sheets{s})(i)-str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1));
                                end
                            end
                        case '-'
                            if j<nsub1 && (alldata{i,1}(psub1(j+1))=='*' || alldata{i,1}(psub1(j+1))=='/')
                                if sum(strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns))
                                    temp=-N([false(1,c1*rxns_l) strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1)))
                                    temp=-str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1));
                                end
                            else
                                if sum(strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns))
                                    L=L-N([false(1,c1*rxns_l) strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1)))
                                    b2.(sheets{s})(i)=b2.(sheets{s})(i)+str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1));
                                end
                            end
                        case '*'
                            if sum(strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns)) && isscalar(temp)
                                L=L+temp*N([false(1,c1*rxns_l) strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns)],:);
                            elseif ~isnan(str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1)))
                                L=L+temp*str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1));
                            else
                                warning([sheets{s} ' ' num2str(i) ': nonlinear constraint'])
                            end
                        case '/'
                            if sum(strcmp(alldata{i,1}(psub1(j)+1:psub1(j+1)-1),rxns))
                                warning([sheets{s} ' ' num2str(i) ': nonlinear constraint'])
                            elseif ~isnan(str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1)))
                                L=L+temp/str2double(alldata{i,1}(psub1(j)+1:psub1(j+1)-1));
                            end
                        otherwise
                            warning([sheets{s} ' ' num2str(i) ': unexpected operation'])
                    end
                end
                for j=1:nsub2
                    switch alldata{i,3}(psub2(j))
                        case '+'
                            if j<nsub2 && (alldata{i,3}(psub2(j+1))=='*' || alldata{i,3}(psub2(j+1))=='/')
                                if sum(strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns))
                                    temp=N([false(1,c1*rxns_l) strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1)))
                                    temp=str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1));
                                end
                            else
                                if sum(strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns))
                                    R=R+N([false(1,c1*rxns_l) strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1)))
                                    b2.(sheets{s})(i)=b2.(sheets{s})(i)-str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1));
                                end
                            end
                        case '-'
                            if j<nsub2 && (alldata{i,3}(psub2(j+1))=='*' || alldata{i,3}(psub2(j+1))=='/')
                                if sum(strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns))
                                    temp=-N([false(1,c1*rxns_l) strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1)))
                                    temp=-str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1));
                                end
                            else
                                if sum(strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns))
                                    R=R-N([false(1,c1*rxns_l) strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns)],:);
                                elseif ~isnan(str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1)))
                                    b2.(sheets{s})(i)=b2.(sheets{s})(i)+str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1));
                                end
                            end
                        case '*'
                            if sum(strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns)) && isscalar(temp)
                                R=R+temp*N([false(1,c1*rxns_l) strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns)],:);
                            elseif ~isnan(str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1)))
                                R=R+temp*str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1));
                            else
                                warning([sheets{s} ' ' num2str(i) ': nonlinear constraint'])
                            end
                        case '/'
                            if sum(strcmp(alldata{i,3}(psub2(j)+1:psub2(j+1)-1),rxns))
                                warning([sheets{s} ' ' num2str(i) ': nonlinear constraint'])
                            elseif ~isnan(str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1)))
                                R=R+temp/str2double(alldata{i,3}(psub2(j)+1:psub2(j+1)-1));
                            end
                        otherwise
                            warning([sheets{s} ' ' num2str(i) ': unexpected operation'])
                    end
                end
                A2.(sheets{s})(i,:)=L-R;
            else
                warning([sheets{s} ' ' num2str(i) ': equality between two numbers'])
            end
        end
    end
    eq.A=[A2.(sheets{1}); A2.(sheets{2})];
    eq.b=[b2.(sheets{1}); b2.(sheets{2})];
    textrep=[textrep2.(sheets{1}); textrep2.(sheets{2})];
    remove_cons=find(sum(eq.A,2)==0);
    for i=1:length(remove_cons)
        if min(eq.A(remove_cons(i),:))==0
            eq.b(remove_cons(i))=0; % remove constraints on nonexistent reactions
        end
    end
end
%% matrix A must have length(net)+length(xch) number of columns
% move constrained portions to the other side
eq.b=eq.b-eq.A(:,rm_cons)*all_free(rm_cons);
eq.A=eq.A(:,~rm_cons);