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

corn_degree_day <- function(x,max_gd,emerg_day){
  for(i in 1:nrow(x)){
    if(x$date[i] < as.Date(emerg_day)){
      x$degree_day[i] <- as.numeric(0.00)
    }else{
      if(x$degree_day[i-1] < max_gd ){
        x$degree_day[i] <- round((x$degree_day[i-1] + (((temp_corn(x$tmax[i]) + temp_corn(x$tmin[i]))/2) - 10)),2)
      }else{
        break
      }
    }
  }
  return(x)
}


test <-bh_2024_2025 %>% 
  select(date,tmax,tmin) %>% 
  corn_degree_day(emerg_day = "2024-10-01",2400)



