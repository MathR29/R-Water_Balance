library(dplyr)

#Segundo universidade de Ohio
# Tmin milho == 10 graus (50F)
# Tmax milho == 30 graus (86F)

temp_corn <- function(number) {
  if (number >= 30) {
    return(30)
  } else if (number >= 10) {
    return(number)
  } else {
    return(10)
  }
}

corn_degree_day <- function(x){
  tmin <- sapply(x$tmin,temp_corn)
  tmax <- sapply(x$tmax,temp_corn)
  x$degree_day <- cumsum(((tmax+tmin)/2)-10)
  return(x)
}

corn_degree_day <- function(x,max_gd,planting_date){
  for(i in 1:nrow(x)){
    if(x$date[i] < as.Date(planting_date)){
      x$degree_day[i] <- as.numeric(0.00)
    }else{
      if(x$degree_day[i-1] < max_gd ){
        x$degree_day[i] <- round((x$degree_day[i-1] + (((temp_corn(x$tmax[i]) + temp_corn(x$tmin[i]))/2) - 10)),2)
      }else{
        break
      }
    }
  }
  x <- x %>% mutate(stage = case_when(degree_day < 100 ~ "-",
                                 degree_day < 346 ~ "VE",
                                 degree_day < 592 ~ "V3",
                                 degree_day < 838 ~ "V6",
                                 degree_day < 1020 ~ "V9",
                                 degree_day < 1170 ~ "V12",
                                 degree_day < 1370 ~ "V15",
                                 degree_day < 1420 ~ "V19",
                                 degree_day < 1686 ~ "VT",
                                 degree_day < 1830 ~ "R2",
                                 degree_day < 1981 ~ "R3",
                                 degree_day < 2324 ~ "R4",
                                 degree_day < 2650 ~ "R5",
                                 TRUE ~ "R6"))
  return(x)
}


test <-bh_2024_2025 %>% 
  select(date,tmax,tmin) %>% 
  corn_degree_day(planting_date = "2024-10-01",2650)



