library(readr)
library(tidyverse)
data <- read.csv("hurrican703.csv")
data <- data %>%
  dplyr::group_by(ID) %>% 
  separate(time, into = c("Date", "Hour"), sep = " ") %>%
  mutate(Date = as.numeric(str_sub(Date,-2,-1)),
         Hour = str_sub(Hour,1,2))%>%
  mutate(Hour = as.numeric(ifelse(Hour == "06", "6", Hour))) %>%
  ungroup() %>%
  mutate(Month = factor(Month,levels = month.name)) %>%
  group_by(ID, Season, Month, Nature, Date) %>%
  mutate(n = n()) %>%
  arrange(Hour, .by_group = TRUE) %>%
  mutate(Wind_change = Wind.kt - lag(Wind.kt, 1),
         Lat_change = Latitude - lag(Latitude, 1),
         Long_change = Longitude - lag(Longitude, 1))
save(data, file = "data.RData")