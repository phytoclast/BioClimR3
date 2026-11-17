
library(climatools)
library(ggplot2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
biosum <- readRDS('data/biomesum.RDS')
    sclim <- biosum[2,]
    sname = paste(sclim$LEVEL_NAME,"|", sclim$probiome)
    mon <- 1:12
    mnn <- sprintf("%02d", mon) 
    lat=sclim$lat
    lon=sclim$lon
    elev=sclim$elev
    df <-  data.frame(mon = 1:12,
                      p=t(sclim[,paste0('p',mnn)]),
                      th=t(sclim[,paste0('th',mnn)]),
                      tl=t(sclim[,paste0('tl',mnn)]))
    
    df$t <- (df$th+df$tl)/2
    
    awc <- 100
    df$Ra <- GetSolarRad(GetDayNumber(df$mon),Lat=lat)
    df$e <- GetPET.df(mon=df$mon,lat=lat,th=df$th,tl=df$tl,p=df$p)
    df<-df|>cbind(SoilWaterBudget(df$p, df$e, awc=awc))
    df$aa <- pmin(df$e, df$p)
    Tg <- climatools::warm6ppt(t=df$t,p=0,n=6,max_var = TRUE, v_input = 't', v_return = 'b',bt=0)
    MAAT <- mean(df$t)
    Tp <- weighted.mean(df$t, df$p)
    Tw <- max(df$t)
    Tc <- min(df$t)
    Tcl <- min(df$tl)
    Twh <- max(df$th)
    PET <- sum(df$e)
    MAP <- sum(df$p)
    AET <- sum(df$aa)
    prat <- MAP/(PET+0.01)
    m <- MAP/(PET+MAP+0.01)*100
    d <- 100-AET/(PET+0.01)*100
    Tclx <- climatools::XtremLow(Tcl=Tcl,Lat=lat,Lon=lon,Elev=elev)
    cindex <- pmin(Tc,Tcl+3)
    cindex2 <- pmin(Tc, Tclx*12/15+12)
    trindex <- Tg-cindex
    climatext <- paste0("Twh: ",round(Twh,1),"°C | Tw: ",round(Tw,1),"°C | Tc: ",round(Tc,1),"°C | Tcl: ",round(Tcl,1),"°C | Tclx: ",round(Tclx,1),"°C\n",
                        "Tg: ",round(Tg,1),"°C | MAAT: ",round(MAAT,1),"°C | Tp: ",round(Tp,1),"°C | PET: ",round(PET,0)," mm\n",
                        "MAP/PET: ",round(prat,2)," | m: ",round(m,0),"% | d: ",round(d,0),"% | MAP: ",round(MAP,0)," mm")
    
    cregime <- case_when(cindex < 0 & Tg - cindex >= 12 ~ 1 ,
                         cindex < 12 & Tg - cindex >= 6  ~ 2,
                         TRUE ~ 3)
    
    mindex <- case_when(prat < 0.125 ~ 'Perarid',
                        prat < 0.25 ~ 'Arid',
                        prat < 0.5 ~ 'Semiarid',
                        prat < 1 ~ 'Subhumid',
                        prat < 2 ~ 'Humid',
                        TRUE ~ 'Perhumid')

    mindex <- case_when(prat < 0.25 ~ mindex,
                        prat < 1 & MAAT > Tp & cregime %in% c(1,2)  ~ paste('Xerothermic',mindex),
                        prat < 1 ~ paste('Pluviothermic',mindex),
                        d >= 25 & MAAT > Tp & cregime %in% c(1,2)  ~ paste('Xerothermic',mindex),
                        d >= 25  ~ paste('Pluviothermic',mindex),
                        TRUE ~ mindex)

    tindex <- case_when(Tg < 6 ~ 'Cold',
                        Tg < 12 ~ 'Cool',
                        Tg < 18 ~ 'Mild',
                        Tg < 24 ~ 'Warm',
                        TRUE ~ 'Hot')

    tindex <- case_when(cregime %in% 1 ~ paste(tindex,'Polar'),
                         cregime %in% 2 ~ paste(tindex,'Temperate'),
                         cregime %in% 3 ~ paste(tindex,'Tropical'))
    
    
    alt <- case_when(cregime  %in% c(1,2) & Tg < 6 ~ 'Arctic/Alpine',
                     cregime  %in% 1 & Tg < 12 ~ 'Boreal/Subalpine',
                     cregime  %in% 1  ~ 'Continental',
                     cregime  %in% 2 & Tg < 18 ~ 'Oceanic',
                     cregime  %in% 2 ~ 'Subtropical',
                     cregime  %in% 3 & Tg < 6 ~ 'Andean',
                     cregime  %in% 3 & Tg < 24  ~ 'Tropical Montane',
                     cregime  %in% 3  ~ 'Tropical')
    
    k1 <- case_when(Tg < 6 ~ 'E',
                    prat < 1 ~ 'B',
                    cregime  %in% 1  ~ 'D',
                    cregime  %in% 2  | Tg < 18 ~ 'C',
                    TRUE  ~ 'A')
    
    
    k2 <- case_when(k1 %in% 'E' & Tw < 0 ~ 'F',
                    k1 %in% 'E'  ~ 'T',
                    k1 %in% 'B' & prat < 0.25  ~ 'W',
                    k1 %in% 'B'  ~ 'S',
                    d >= 25 & MAAT > Tp  ~ 's',
                    d >= 25   ~ 'w',
                    TRUE  ~ 'f')
    k3 <- case_when(k1 %in% c('E','A') ~ '',
                    k1 %in% 'B' & Tg < 18   ~ 'k',
                    k1 %in% 'B'  ~ 'h',
                    k1 %in% c('C','D') & Tg < 12 ~ 'c',
                    k1 %in% c('C','D') & Tg < 18 ~ 'b',
                    TRUE ~ 'a')


    climaclass <- paste(sname,'\n',mindex, tindex, paste0(' ',k1,k2,k3, ' (',alt,')'))
    
    
    # rv$climatext <- climatext

    makeClimplot( t=df$t,
                  p=df$p,
                  th = df$th,
                  tl = df$tl,
                  e = df$e,
                  a = df$a,
                  u = df$u,
                  t80 = NULL,
                  t20 = NULL,
                  p80 = NULL,
                  p20 = NULL,
                  name = climaclass,
                  bar = F)
    
    a1=data.frame(x=c(-50,-50,0,-6), y=c(0,6,6,0))
    a2=data.frame(x=c(-50,-50,0,-6), y=c(6,12,12,6))
    a3=data.frame(x=c(-50,-50,0,0), y=c(12,36,36,12))
    a4=data.frame(x=c(-6,0,6,0), y=c(0,6,6,0))
    a5=data.frame(x=c(-6,0,0,12,0), y=c(6,12,18,18,6))
    a6=data.frame(x=c(0,0,12,12), y=c(18,36,36,18))
    a7=data.frame(x=c(0,12,12,24,18,6), y=c(6,18,24,24,18,6))
    a8=data.frame(x=c(12,12,36,24), y=c(24,36,36,24))
    breaksx=c(-45,-40, -35, -30, -25, -20,-15, -10,-5, 0,5, 10,15, 20,25,30)
    breaksy=c(0,6,12,18,24,30)
    ll1 <- data.frame(x=c(-50,6), y=c(6,6))
    ll2 <- data.frame(x=c(0,0), y=c(0,36))
    ll3 <- data.frame(x=c(12,12), y=c(18,36))
    l1 <- data.frame(x=c(-50,12), y=c(12,12))
    l2 <- data.frame(x=c(-12,0), y=c(0,12))
    l3 <- data.frame(x=c(-50,18), y=c(18,18))
    l4 <- data.frame(x=c(-50,24), y=c(24,24))
    l5 <- data.frame(x=c(-6,12), y=c(0,18))
    l6 <- data.frame(x=c(-25,-25), y=c(15,36))
    l7 <- data.frame(x=c(5,5), y=c(18,36))
    l8 <- data.frame(x=c(12,12), y=c(24,36))
    # climplot2 <-  
    ggplot() +
      geom_polygon(data=a1, mapping=aes(x=x, y=y, fill='alpine'),alpha = 0.5)+
      geom_polygon(data=a2, mapping=aes(x=x, y=y, fill='boreal'),alpha = 0.5)+
      geom_polygon(data=a3, mapping=aes(x=x, y=y, fill='continental'),alpha = 0.5)+
      geom_polygon(data=a4, mapping=aes(x=x, y=y, fill='andean'),alpha = 0.5)+
      geom_polygon(data=a5, mapping=aes(x=x, y=y, fill='oceanic'),alpha = 0.5)+
      geom_polygon(data=a6, mapping=aes(x=x, y=y, fill='subtropical'),alpha = 0.5)+
      geom_polygon(data=a7, mapping=aes(x=x, y=y, fill='montane'),alpha = 0.5)+
      geom_polygon(data=a8, mapping=aes(x=x, y=y, fill='tropical'),alpha = 0.5)+
      
      geom_line(data=ll1, mapping=aes(x=x, y=y),alpha = 0.2, color='black', linetype='solid')+
      geom_line(data=ll2, mapping=aes(x=x, y=y),alpha = 0.2, color='black', linetype='solid')+
      geom_line(data=ll3, mapping=aes(x=x, y=y),alpha = 0.2, color='black', linetype='solid')+
      geom_line(data=l1, mapping=aes(x=x, y=y),alpha = 0.2, color='black', linetype='solid')+
      geom_line(data=l2, mapping=aes(x=x, y=y),alpha = 0.2, color='black', linetype='solid')+
      geom_line(data=l3, mapping=aes(x=x, y=y),alpha = 0.2, color='black', linetype='solid')+
      geom_line(data=l4, mapping=aes(x=x, y=y),alpha = 0.2, color='black', linetype='solid')+
      geom_line(data=l5, mapping=aes(x=x, y=y),alpha = 0.2, color='black', linetype='solid')+
      # geom_line(data=l6, mapping=aes(x=x, y=y),alpha = 0.2, color='black', linetype='solid')+
      # geom_line(data=l7, mapping=aes(x=x, y=y),alpha = 0.2, color='black', linetype='solid')+
      geom_line(data=l8, mapping=aes(x=x, y=y),alpha = 0.2, color='black', linetype='solid')+
      # geom_point(data=sclim, mapping=aes(x=cindex, y=Tg), color = 'black', size=0.5)+
      #geom_density2d(data=StationMeans, mapping=aes(x=Cindex, y=Tg),color = 'black',alpha = 0.25)+
      scale_fill_manual("Legend", values = c("alpine" = "lightblue",
                                             "boreal" = "darkgreen",
                                             "continental" = "yellowgreen",
                                             "andean" = "yellow",
                                             "oceanic" = "darkcyan",
                                             "subtropical" = "orange",
                                             "montane" = "pink",
                                             "tropical" = "darkred"
      ))+
      
      scale_x_continuous(name= "Coldest Month (Annual Extreme Minimum)",
                         breaks=breaksx,
                         labels=paste0(breaksx,'(', breaksx-3,')'))+
      scale_y_continuous(name= "Growing Season", breaks=breaksy, expand=F)+
      coord_fixed(ratio = 1/1,xlim = c(-45,30), ylim = c(0, 33))+
      labs(title = paste("Climate of ",sname, "", sep=""))+
      theme_bw()+
      theme(legend.position='right',axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0),
            panel.grid.major = element_line(), panel.grid.minor = element_blank())