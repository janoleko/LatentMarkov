library(tidyverse)
setwd("/Users/jan-ole/R/Packages_on_git/LatentMarkov")

all_handball_files = list.files(path = "./data/handball_matches_2022-23/", pattern = ".csv")
all_matches = paste0("./data/handball_matches_2022-23/", all_handball_files) %>% 
  map(~read.csv(., header = TRUE))

# creating a matchID
for(i in 1:length(all_matches)){ 
  all_matches[[i]]$matchID = i
  if(is.integer(all_matches[[i]][,1])){
    all_matches[[i]] = all_matches[[i]][,-1]
  }
}

head(all_matches[[2]])

# some cleanup
all_matches[[97]]$Penalty_Hometeam[which(all_matches[[97]]$Penalty_Hometeam=="o")] = "0"
all_matches[[97]]$Penalty_Hometeam = all_matches[[97]]$Penalty_Hometeam %>% as.integer()

# binding all matches to one data frame
all_matches = bind_rows(all_matches)

# omit games without time variable
all_matches = all_matches %>% filter(!is.na(Time))

# more cleaning (there are two columns for the team names)
all_matches$hometeam[which(is.na(all_matches$hometeam))] = all_matches$Hometeam[which(is.na(all_matches$hometeam))]
all_matches$awayteam[which(is.na(all_matches$awayteam))] = all_matches$Awayteam[which(is.na(all_matches$awayteam))]

all_matches$Hometeam[which(is.na(all_matches$Hometeam))] = all_matches$hometeam[which(is.na(all_matches$Hometeam))]
all_matches$Awayteam[which(is.na(all_matches$Awayteam))] = all_matches$awayteam[which(is.na(all_matches$Awayteam))]

# deselecting the columns that are not needed
all_matches = all_matches %>% select(-Hometeam, -Awayteam)

# fill in blank scores
all_matches$Score[which(all_matches$Score=="")] = NA
all_matches = all_matches %>% fill(Score, .direction = "up")

# create home and away goals variables
all_matches$home_goals = as.integer(str_extract(all_matches$Score, "^[0-9]+"))
all_matches$away_goals = as.integer(str_extract(all_matches$Score, "[0-9]+$"))


# Find 7-meter strings ----------------------------------------------------

str_detect(all_matches$Event.Description, "7-Meter") %>% sum()
# 2064 7-meter throws

# Find 7-meter attempts and create success variable
all_matches$seven_meter_attempt = as.integer(str_detect(all_matches$Event.Description, "7-Meter"))

attemptInd = which(all_matches$seven_meter_attempt==1)

succesInd = attemptInd[which(all_matches$home_goals[attemptInd] != all_matches$home_goals[attemptInd+1] | 
                               all_matches$away_goals[attemptInd] != all_matches$away_goals[attemptInd+1])]

all_matches$seven_meter_success = 0
all_matches$seven_meter_success[succesInd] = 1



# Filter out everything else ----------------------------------------------

data = all_matches %>% filter(seven_meter_attempt==1) %>% select(matchID, time = Time,
                                                                 hometeam, awayteam, 
                                                                 home_goals, away_goals, 
                                                                 seven_meter_success, event = Event.Description)

# first detect brackets in Event
brackets = which(str_detect(data$event, "\\(.*\\)"))
data$event[brackets]
# extract player names (two words) in front of brackets
data$player = NA
data$player[brackets] = str_extract(data$event[brackets], "\\w+\\s\\w+(?= \\()")

NAind = which(is.na(data$player))
unique_players = unique(data$player[brackets])
# find the player names in players vector in the remaining events
for(i in 1:length(unique_players)){
  data$player[which(str_detect(data$event, unique_players[i]))] = unique_players[i]
}

NAind = which(is.na(data$player))
data$event[NAind]
data$player[NAind][3] = "Julian Köster"

# formatting time variable
str_length(data$time)
longInd = which(str_length(data$time)>5)
data$time[longInd] = str_sub(data$time[longInd], 1, 5)
data$time2 = strptime(data$time, format = "%M:%S")

# reordering columns
data = data %>% select(matchID, time, time2, hometeam, awayteam, home_goals, away_goals, seven_meter_success, player, event)

# reordering for unique players
unique_players = unique(data$player)
data2 = list()
for(i in 1:length(unique_players)){
  data2[[i]] = data %>% filter(player == unique_players[i])
}

data2 = bind_rows(data2)

# save data

write.csv(data2, file = "./data/handball_data_7m.csv", row.names = FALSE)

