## Livrarias / Pacotes ##
library(tidyverse)
library(nasapower)


## Criando função para obtenção dos dados climáticos
dados_climaticos_cidades <- function(path,data_inicio,data_final){
  df_cidades <- read.csv(path)
  dados_climaticos <- data.frame()
  for (i in 1:nrow(df_cidades)){
    df <- get_power(
      pars = c("T2M","T2M_MAX","T2M_MIN","PRECTOTCORR"),
      community = "AG",
      temporal_api = "daily",
      lonlat = c(df_cidades$long[i],df_cidades$lat[i]),
      dates = c(data_inicio,data_final)) %>% 
      rename(
        date = YYYYMMDD,
        latitude = LAT,
        longitude = LON,
        tmean = T2M, # TEMPERATURA MÉDIA
        tmax = T2M_MAX, # TEMPERATURA MÁXIMA
        tmin = T2M_MIN, # TEMPERATURA MÍNIMA
        rain = PRECTOTCORR, # PRECIPITAÇÃO
      ) %>% 
      mutate(cidade = df_cidades$cidades[i],
             .before = longitude) %>% 
      mutate(
        cad = df_cidades$cad[i],
        lat_rad = latitude*0.0174533,
        dr = 1 + 0.033*cos((2*pi/365)*DOY),
        Sd = 0.409*sin((2*pi/365)*DOY - 1.39),
        ws = acos(-tan(lat_rad)*tan(Sd)),
        Ra = (24*60)/(pi) * 0.0820 * dr * (ws*sin(lat_rad)*sin(Sd)+ cos(lat_rad)*sin(ws)),
        etp = 0.0135 * 0.17  * (Ra / 2.45) * (sqrt(tmax-tmin)) * (tmean + 17.8),
        p_etp = (rain - etp)) %>% 
      select(-c(dr, Sd, ws, Ra,lat_rad)) %>% 
      group_by(cidade) %>% 
      nest() %>% 
      ungroup()
    
    dados_climaticos <- bind_rows(dados_climaticos,df)
  }
  return(dados_climaticos)
}

## Função para cálculo do balanço hídrico
balanco_hidrico <- function(path,data_inicio,data_final){
  dados_climaticos <- dados_climaticos_cidades(path,data_inicio,data_final)
  for (i in 1:nrow(dados_climaticos)){
    for (j in 1:nrow(dados_climaticos$data[[i]])){
      if (j == 1){
        if (dados_climaticos$data[[i]]$p_etp[j] < 0){
          dados_climaticos$data[[i]]$nac[j] <- (dados_climaticos$data[[i]]$p_etp[j])
          dados_climaticos$data[[i]]$arm[j] <- dados_climaticos$data[[i]]$cad[j] * exp((dados_climaticos$data[[i]]$nac[j]) / dados_climaticos$data[[i]]$cad[j])
        }
        if (dados_climaticos$data[[i]]$p_etp[j] >= 0){
          dados_climaticos$data[[i]]$arm[j] <- dados_climaticos$data[[i]]$cad[j] + dados_climaticos$data[[i]]$p_etp[j]
          dados_climaticos$data[[i]]$nac[j] <- dados_climaticos$data[[i]]$cad[j] * log(dados_climaticos$data[[i]]$arm[j] / dados_climaticos$data[[i]]$cad[j])
        }
      }else{
        if (dados_climaticos$data[[i]]$p_etp[j] < 0){
          dados_climaticos$data[[i]]$nac[j] <- (dados_climaticos$data[[i]]$p_etp[j] + dados_climaticos$data[[i]]$nac[j-1]  )
          dados_climaticos$data[[i]]$arm[j] <- dados_climaticos$data[[i]]$cad[j] * exp((dados_climaticos$data[[i]]$nac[j]) / dados_climaticos$data[[i]]$cad[j])
        }
        if (dados_climaticos$data[[i]]$p_etp[j] >= 0){
          dados_climaticos$data[[i]]$arm[j] <- dados_climaticos$data[[i]]$arm[j-1]  + dados_climaticos$data[[i]]$p_etp[j]
          dados_climaticos$data[[i]]$nac[j] <- dados_climaticos$data[[i]]$cad[j] * log(dados_climaticos$data[[i]]$arm[j] / dados_climaticos$data[[i]]$cad[j])
        }
      }
      if (dados_climaticos$data[[i]]$nac[j] >= 0){
        dados_climaticos$data[[i]]$nac[j] <- 0
      }
      if (dados_climaticos$data[[i]]$arm[j] > dados_climaticos$data[[i]]$cad[j]){
        dados_climaticos$data[[i]]$arm[j] <- dados_climaticos$data[[i]]$cad[j]
      }
    }
    dados_climaticos$data[[i]] <- mutate(.data = dados_climaticos$data[[i]],
                                         alt =c(0, diff(arm)),
                                         etr = case_when(p_etp < 0 ~ rain + abs(alt),
                                                         p_etp >= 0 ~ etp),
                                         def = etp - etr,
                                         exc =case_when(arm < cad ~ 0,
                                                        arm == cad ~ p_etp - alt))
  }
  return(dados_climaticos %>% unnest(data))
}

## Função para cálculo do balanço hídrico da cultura da soja
balanco_hidrico_soja <- function(path,data_inicio,data_final,data_emergencia,ciclo){
  df <- balanco_hidrico(path,data_inicio,data_final) %>% 
    group_by(cidade) %>% 
    nest() %>% 
    ungroup()
  for (i in 1:nrow(df)){
    for (j in 1:nrow(df$data[[i]])){
      df$data[[i]]$DAE[j] <- case_when(
        df$data[[i]]$date[j] < as.Date(data_emergencia) ~ NA,
        df$data[[i]]$date[j] >= as.Date(data_emergencia) & df$data[[i]]$date[j] < (as.Date(data_emergencia) + ciclo) ~  
          as.numeric(df$data[[i]]$date[j] - as.Date(data_emergencia)) + 1)
      
    }
    if (ciclo == 110){
      df$data[[i]] <-  df$data[[i]] %>% inner_join(kc_values_110) %>% 
        mutate(etc = etp * kc,
               p_etc = rain - etc) %>% 
        select(-c(arm,nac,def,etr,exc,alt,p_etp))
    }
    if (ciclo == 120){
      df$data[[i]] <-  df$data[[i]] %>% inner_join(kc_values_120) %>% 
        mutate(etc = etp * kc,
               p_etc = rain - etc) %>% 
        select(-c(arm,nac,def,etr,exc,alt,p_etp))
    }
    if (ciclo == 130){
      df$data[[i]] <-  df$data[[i]] %>% inner_join(kc_values_130) %>% 
        mutate(etc = etp * kc,
               p_etc = rain - etc) %>% 
        select(-c(arm,nac,def,etr,exc,alt,p_etp))
    }
  }
  for (i in 1:nrow(df)){
    for (j in 1:nrow(df$data[[i]])){
      if (j == 1){
        if (df$data[[i]]$p_etc[j] < 0){
          df$data[[i]]$nac[j] <- (df$data[[i]]$p_etc[j])
          df$data[[i]]$arm[j] <- df$data[[i]]$cad[j] * exp((df$data[[i]]$nac[j]) / df$data[[i]]$cad[j])
        }
        if (df$data[[i]]$p_etc[j] >= 0){
          df$data[[i]]$arm[j] <- df$data[[i]]$cad[j] + df$data[[i]]$p_etc[j]
          df$data[[i]]$nac[j] <- df$data[[i]]$cad[j] * log(df$data[[i]]$arm[j] / df$data[[i]]$cad[j])
        }
      }else{
        if (df$data[[i]]$p_etc[j] < 0){
          df$data[[i]]$nac[j] <- (df$data[[i]]$p_etc[j] + df$data[[i]]$nac[j-1]  )
          df$data[[i]]$arm[j] <- df$data[[i]]$cad[j] * exp((df$data[[i]]$nac[j]) / df$data[[i]]$cad[j])
        }
        if (df$data[[i]]$p_etc[j] >= 0){
          df$data[[i]]$arm[j] <- df$data[[i]]$arm[j-1]  + df$data[[i]]$p_etc[j]
          df$data[[i]]$nac[j] <- df$data[[i]]$cad[j] * log(df$data[[i]]$arm[j] / df$data[[i]]$cad[j])
        }
      }
      if (df$data[[i]]$nac[j] >= 0){
        df$data[[i]]$nac[j] <- 0
      }
      if (df$data[[i]]$arm[j] > df$data[[i]]$cad[j]){
        df$data[[i]]$arm[j] <- df$data[[i]]$cad[j]
      }}
    df$data[[i]] <- mutate(.data = df$data[[i]],
                           alt =c(0, diff(arm)),
                           etr = case_when(p_etc < 0 ~ rain + abs(alt),
                                           p_etc >= 0 ~ etc),
                           def = etc - etr,
                           exc =case_when(arm < cad ~ 0,
                                          arm == cad ~ p_etc - alt),
                           `etr/etc` = etr/etc)
  }
  return(df %>% unnest(data))
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
    facet_wrap(~cidade)+
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
      facet_wrap(~cidade)+
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
      facet_wrap(~cidade)+
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
    facet_wrap(~cidade)+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  return(plot)
}




