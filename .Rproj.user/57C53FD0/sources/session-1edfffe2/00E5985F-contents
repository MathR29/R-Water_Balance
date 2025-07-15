## Livrarias / Pacotes ##
library(tidyverse)
library(nasapower)


## Criação da função para calculo do BH diário ##
balanco_hidrico <- function(lat,lon,data_inicio,data_final,cad){ 
  dados_climaticos <- get_power(
    pars = c("T2M","T2M_MAX","T2M_MIN","PRECTOTCORR"),
    community = "AG",
    temporal_api = "daily",
    lonlat = c(lon,lat),
    dates = c(data_inicio,data_final))  %>% 
    rename(
      date = YYYYMMDD,
      latitude = LAT,
      longitude = LON,
      tmean = T2M, # TEMPERATURA MÉDIA
      tmax = T2M_MAX, # TEMPERATURA MÁXIMA
      tmin = T2M_MIN, # TEMPERATURA MÍNIMA
      rain = PRECTOTCORR, # PRECIPITAÇÃO
    ) %>% 
    mutate(
      lat_rad = latitude*0.0174533,
      dr = 1 + 0.033*cos((2*pi/365)*DOY),
      Sd = 0.409*sin((2*pi/365)*DOY - 1.39),
      ws = acos(-tan(lat_rad)*tan(Sd)),
      Ra = (24*60)/(pi) * 0.0820 * dr * (ws*sin(lat_rad)*sin(Sd)+ cos(lat_rad)*sin(ws)),
      etp = 0.0135 * 0.17  * (Ra / 2.45) * (sqrt(tmax-tmin)) * (tmean + 17.8),
      p_etp = (rain - etp)) %>% 
    select(-c(dr, Sd, ws, Ra,lat_rad))
  
  
  
  for (i in 1:nrow(dados_climaticos)) { ##For loop para calcular a ETP 
    if (i == 1){
      if (dados_climaticos$p_etp[i] < 0) {
        
        dados_climaticos$nac[i] <- (dados_climaticos$p_etp[i])
        dados_climaticos$arm[i] <- cad * exp((dados_climaticos$nac[i]) / cad)
      }
      if (dados_climaticos$p_etp[i] >= 0){
        dados_climaticos$arm[i] <- cad + dados_climaticos$p_etp[i]
        dados_climaticos$nac[i] <- cad * log(dados_climaticos$arm[i] / cad)
      }
    }
    if (i != 1){
      if (dados_climaticos$p_etp[i] < 0) {
        
        dados_climaticos$nac[i] <- (dados_climaticos$p_etp[i] + dados_climaticos$nac[i-1])
        
        dados_climaticos$arm[i] <- cad * exp((dados_climaticos$nac[i]) / cad)
      }
      if (dados_climaticos$p_etp[i] >= 0){
        dados_climaticos$arm[i] <- dados_climaticos$arm[i-1] + dados_climaticos$p_etp[i]
        dados_climaticos$nac[i] <- cad * log(dados_climaticos$arm[i] / cad)
      }
    }
    if (dados_climaticos$nac[i] >= 0){
      dados_climaticos$nac[i] <- 0
    }
    if (dados_climaticos$arm[i] > cad){
      dados_climaticos$arm[i] <- cad
    }
  }
  
  
  
  return(dados_climaticos %>%   
           mutate(cad = cad,
                  alt = c(0, diff(arm)),
                  etr = case_when(p_etp < 0 ~ rain + abs(alt),
                                  p_etp >= 0 ~ etp),
                  def = etp - etr,
                  exc = case_when(arm < cad ~ 0,
                                  arm == cad ~ p_etp - alt)
           )
  )
}

## #Criação da função para calculo do BH diário para cultura da soja ##
balanco_hidrico_soja <- function(df_balanco_hidrico,dia_da_emergencia,ciclo){
  if (ciclo == 110){
    for (i in 1:nrow(df_balanco_hidrico)) { 
      
      df_balanco_hidrico$DAE[i] <- case_when(
        df_balanco_hidrico$date[i] < as.Date(dia_da_emergencia) ~ NA,
        df_balanco_hidrico$date[i] >= as.Date(dia_da_emergencia) & df_balanco_hidrico$date[i] < (as.Date(dia_da_emergencia)+110) ~ 
          as.numeric(df_balanco_hidrico$date[i] - as.Date(dia_da_emergencia)) + 1)
      
    }
    df_balanco_hidrico <- df_balanco_hidrico %>% 
      inner_join(kc_values_110) 
  }
  
  if (ciclo == 120){
    for (i in 1:nrow(df_balanco_hidrico)) { 
      
      df_balanco_hidrico$DAE[i] <- case_when(
        df_balanco_hidrico$date[i] < as.Date(dia_da_emergencia) ~ NA,
        df_balanco_hidrico$date[i] >= as.Date(dia_da_emergencia) & df_balanco_hidrico$date[i] < (as.Date(dia_da_emergencia)+120) ~ 
          as.numeric(df_balanco_hidrico$date[i] - as.Date(dia_da_emergencia)) + 1)
      
    }
    df_balanco_hidrico <- df_balanco_hidrico %>% 
      inner_join(kc_values_120) 
  }
  
  if (ciclo == 130){
    for (i in 1:nrow(df_balanco_hidrico)) { 
      
      df_balanco_hidrico$DAE[i] <- case_when(
        df_balanco_hidrico$date[i] < as.Date(dia_da_emergencia) ~ NA,
        df_balanco_hidrico$date[i] >= as.Date(dia_da_emergencia) & df_balanco_hidrico$date[i] < (as.Date(dia_da_emergencia)+130) ~ 
          as.numeric(df_balanco_hidrico$date[i] - as.Date(dia_da_emergencia)) + 1)
      
    }
    df_balanco_hidrico <- df_balanco_hidrico %>% 
      inner_join(kc_values_130) 
  }
  df_balanco_hidrico <- df_balanco_hidrico%>% 
    mutate(etc = etp * kc,
           p_etc = rain - etc) %>% 
    select(-c(arm,nac,def,etr,exc,alt,p_etp))
  
  
  for (i in 1:nrow(df_balanco_hidrico)) { ##For loop para calcular a ETC
    if (i == 1){
      if (df_balanco_hidrico$p_etc[i] < 0) {
        
        df_balanco_hidrico$nac[i] <- (df_balanco_hidrico$p_etc[i])
        df_balanco_hidrico$arm[i] <- df_balanco_hidrico$cad[i] * exp((df_balanco_hidrico$nac[i]) / df_balanco_hidrico$cad[i])
      }
      if (df_balanco_hidrico$p_etc[i] >= 0){
        df_balanco_hidrico$arm[i] <- df_balanco_hidrico$cad[i] + df_balanco_hidrico$p_etc[i]
        df_balanco_hidrico$nac[i] <- df_balanco_hidrico$cad[i] * log(df_balanco_hidrico$arm[i] / df_balanco_hidrico$cad[i])
      }
    }
    if (i != 1){
      if (df_balanco_hidrico$p_etc[i] < 0) {
        
        df_balanco_hidrico$nac[i] <- (df_balanco_hidrico$p_etc[i] + df_balanco_hidrico$nac[i-1])
        
        df_balanco_hidrico$arm[i] <- df_balanco_hidrico$cad[i] * exp((df_balanco_hidrico$nac[i]) / df_balanco_hidrico$cad[i])
      }
      if (df_balanco_hidrico$p_etc[i] >= 0){
        df_balanco_hidrico$arm[i] <- df_balanco_hidrico$arm[i-1] + df_balanco_hidrico$p_etc[i]
        df_balanco_hidrico$nac[i] <- df_balanco_hidrico$cad[i] * log(df_balanco_hidrico$arm[i] / df_balanco_hidrico$cad[i])
      }
    }
    if (df_balanco_hidrico$nac[i] >= 0){
      df_balanco_hidrico$nac[i] <- 0
    }
    if (df_balanco_hidrico$arm[i] > df_balanco_hidrico$cad[i]){
      df_balanco_hidrico$arm[i] <- df_balanco_hidrico$cad[i]
    }
  }
  
  return(df_balanco_hidrico %>% 
           mutate(cad = df_balanco_hidrico$cad[i],
                  alt = c(0, diff(arm)),
                  etr = case_when(p_etc < 0 ~ rain + abs(alt),
                                  p_etc >= 0 ~ etc),
                  def = etc - etr,
                  exc = case_when(arm < df_balanco_hidrico$cad[i] ~ 0,
                                  arm == df_balanco_hidrico$cad[i] ~ p_etc - alt),
                  `etr/etc` = etr/etc))
}

## Criando funçoes para gerar os graficos ##
# Excesso e defict
g_exc_def <- function(df){
  df_long <- df %>% 
    mutate(def = def * -1) %>% 
    pivot_longer(cols = c(exc, def), 
                 names_to = "Tipo", 
                 values_to = "Chuva")
  
  plot <- ggplot(df_long, aes(x = date, y = Chuva, color = Tipo, group = Tipo)) +
    geom_line(size = 1) +
    scale_x_date(date_breaks = "1 week", date_labels = "%d/%m") +
    labs(x = "Mês", y = "mm", title = "Extrato do Balaço Hídrico Mensal Normal") +
    scale_color_manual(name = "", values = c("exc" = "blue", "def" = "red")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  return(plot)
}

# Chuva, Evapotranspiracao potencial ou da cultura e Evapotranspiracao real
g_rain_etp_etr <- function(df,tipo_bh){
  if (tipo_bh == "sequencial") {
    df_long <- df  %>%
      pivot_longer(cols = c(rain,etp, etr), 
                   names_to = "Tipo", 
                   values_to = "Chuva")
    
    plot <- ggplot(df_long, aes(x = date, y = Chuva, color = Tipo, group = Tipo)) +
      geom_line(size = 1) +
      labs(x = "Mês", y = "mm", title = "Balaço Hídrico Normal Mensal") +
      scale_x_date(date_breaks = "1 week", date_labels = "%d/%m")+
      scale_color_manual(name = "", values = c("rain" = "blue", "etp" = "red","etr" = "green" )) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5))
    return(plot)
  }
  if (tipo_bh == "cultura"){
    df_long <- df  %>%
      pivot_longer(cols = c(rain,etc, etr), 
                   names_to = "Tipo", 
                   values_to = "Chuva")
    
    plot <- ggplot(df_long, aes(x = date, y = Chuva, color = Tipo, group = Tipo)) +
      geom_line(size = 1) +
      scale_x_date(date_breaks = "1 week", date_labels = "%d/%m")+
      labs(x = "Mês", y = "mm", title = "Balaço Hídrico Normal Mensal") +
      scale_color_manual(name = "", values = c("rain" = "blue", "etc" = "red","etr" = "green" )) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5))
    return(plot)
  }
}

# Capacidade de armazenamento de agua e Armazenamento
g_cad_arm <- function(df){
  df_long <- df %>%
    pivot_longer(cols = c(cad, arm), 
                 names_to = "Tipo", 
                 values_to = "Chuva")
  plot <-  ggplot(df_long, aes(x = date, y = Chuva, color = Tipo, group = Tipo)) +
    geom_line(size = 1) +
    labs(x = "Mês", y = "mm", title = "Capacidade de Armazenamento (CAD) e Armazenamento (Arm) mensal ") +
    scale_x_date(date_breaks = "1 week", date_labels = "%d/%m") +
    scale_color_manual(name = "", values = c("cad" = "blue", "arm" = "red")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  return(plot)
}




