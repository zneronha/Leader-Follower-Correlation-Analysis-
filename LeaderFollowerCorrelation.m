function LeaderFollowerCorrelation
clearvars; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Version 1.0 2 January 2018: basic functionality, bug fixes
%Version 2.0 5 January 2018: computational completion
%Version 2.1 5 January 2018: data storage improvements
%Version 3.0 9 January 2018: multi-leadeer follower searching bug fixed 
%Version 3.1 11 January 2018: detailed commenting, saving feature

%load compiled dataset 
load('Z:\ENG_BBCancer_Shared\group\0Zach\LeaderCombinedData\EGF(E6)_compiledleaderdata.mat',...
    'CompiledDataStore');
%extract relevant well data from the array
u = ~cellfun('isempty',CompiledDataStore); %#ok<USENS>
[wells,~,~] = find(u(:,1)==1);
frames = 1:12;

%Allocate memory for total storage cell
MacroStore = cell(max(wells),max(frames));

%initialize wellcounter
wellcounter = 0;
for well = wells'
    %sequence WC and initialize frame counter
    wellcounter = wellcounter + 1;
    framecounter = 0;
    %store well data for particular well in a temporary matrix for
    %referencing
    tempwelldata = CompiledDataStore{well,1}; %#ok<IDISVAR>

    %load velocity data
    loadname = strcat('Z:\ENG_BBCancer_Shared\group\0Zach\Cluster data\EGF(E3E6)LeaderVelData\vel-EGF(E6)w',...
        num2str(well),'.mat');
    load(loadname,'velX','velY');
    
    %loop through the frames
    for frame = frames
        %sequence frame counter and initilize usercounter
        framecounter = framecounter + 1;
        personcounter = 0;
        
        %reserve memory for data storage matrix
        MacroUserStore = cell(1,3);
        
        %loop through each one of the users
        for user = 1:3
            personcounter = personcounter + 1;
            
            %store leader and follower data for user
            tempFPLeaders = tempwelldata{2*user-1,frame};
            tempFPFollowers = tempwelldata{2*user,frame};
            
            %compute number of leaders
            numleaders = sum(~isnan(tempFPLeaders(:,1)));
            
            %safety check for null frames, prevents code from crashing
            if numleaders == 0
                continue
            end
            
            %allocate memory for this user's data storage
            tempdatastoreLA = cell(numleaders,3);
            
            for leader = 1:numleaders

                %load the leader ID
                leaderID = tempFPLeaders(leader,1);
                tvx = velX(leaderID,frame*5);
                tvy = velY(leaderID,frame*5);
                LeaderUnitV = [tvx tvy]./sqrt(tvx^2+tvy^2);
                
                
                %HERE WE NEED TO IDENTIFY THE CORRECT FOLLOWERS SET
                for qualcheck = 1:numleaders %run through the possible follower sets
                    
                    %load and quantify data size
                   qualfollowermat = tempFPFollowers{qualcheck,1}; 
                   numfollowers = sum(~isnan(qualfollowermat(:,1)));
                   
                   %throw an error if no followers are identified
                    if numfollowers<1
                       errormessage = strcat({'CRASH: A LEADER HAS BEEN IDENTIFIED BY USER'},...
                           num2str(user),{'IN WELL'},num2str(well),{'FRAME'},num2str(frame));
                       error(errormessage)
                    end
                
                    %loop through the followers to determine if there is an
                    %app. match
                   for follower = 1:numfollowers
                       followerID = qualfollowermat(follower,1);
                       
                       %if so, YAY!
                       if followerID == leaderID 
                          QUALFOLLOWER = qualcheck;
                           continue
                       end   
                       %if not, continue...
                   end
                    
                end
                
                followersMAT = tempFPFollowers{QUALFOLLOWER,1};
                numfollowers = sum(~isnan(followersMAT(:,1)));
                
                %throw an error if no followers are identified
                if numfollowers<1
                   errormessage = strcat({'CRASH: A LEADER HAS BEEN IDENTIFIED BY USER'},...
                       num2str(user),{'IN WELL'},num2str(well),{'FRAME'},num2str(frame),...
                       {'WITH NO FOLLOWERS!'});
                   error(errormessage)
                end
                
                %allocate memory for data storage
                tempfollowerstore = nan(numfollowers,2);
                
                %loop through the followers
                for follower = 1:numfollowers
                    %pull the follower ID on this particular follower
                    followerID = followersMAT(follower,1);
                    
                    %Verify comparisons are only made for distinct cells in
                    %the complex
                   
                    if followerID == leaderID
                        continue
                    end
                    %pull velocities, find unit vector
                    fvx = velX(followerID,frame*5);
                    fvy = velY(followerID,frame*5);
                    FollowerUnitV = [fvx fvy]./sqrt(fvx^2+fvy^2);
                    
                    %calculate correlation coefficient and store it with
                    %the cell ID
                    Correlation = dot(LeaderUnitV,FollowerUnitV);
                    tempfollowerstore(follower,:) = [followerID Correlation];
                end
                
                %calculate means and store in appropriate location
                meanCorr = nanmean(tempfollowerstore(:,2));
                
                %store the data in the appropriate location
                tempdatastoreLA{leader,1} = leaderID;
                tempdatastoreLA{leader,2} = tempfollowerstore;
                tempdatastoreLA{leader,3} = meanCorr;
                
            end
            %store the leader cell data
            MacroUserStore{1,user} = tempdatastoreLA;
        end
        %store the well data
        MacroStore{well,frame} = MacroUserStore;
    end
end

save('L-F_Correlation.mat','MacroStore');
end