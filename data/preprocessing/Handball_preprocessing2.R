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

# creating score_home variable
all_matches$score_home = as.integer(all_matches$Score %>% str_extract("^[0-9]+"))
# creating score away variable
all_matches$score_away = as.integer(all_matches$Score %>% str_extract("[0-9]+$"))

# filling in blanks
all_matches = all_matches %>% fill(score_home, .direction = "up")
all_matches = all_matches %>% fill(score_away, .direction = "up")

# creating variable if home team scored
all_matches$home_scored = 0
all_matches$home_scored[which(diff(all_matches$score_home) == -1)] = 1

# creating variable if away team scored
all_matches$away_scored = 0
all_matches$away_scored[which(diff(all_matches$score_away) == -1)] = 1

# reordering columns
handball_data = all_matches %>% 
  select(matchID, Time, hometeam, awayteam, score_home, score_away, home_scored, away_scored, event = Event.Description)

# reverse order
n = nrow(handball_data)
handball_data = handball_data[n:1,]

# save data
write.csv(handball_data, "./data/handball_data.csv", row.names = FALSE)
