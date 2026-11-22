
library(climatools)
library(ggplot2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

biomepts <- readRDS('data/biomepts.RDS')
summits <- readRDS('data/summits.RDS')
summits2 <- summits
biomepts <- biomepts |> rbind(summits[,colnames(biomepts)])
mon <- 1:12
mnn <- sprintf("%02d", mon) 
thesecols <- c(paste0('th',mnn),paste0('th',mnn))
for(i in 1:nrow(summits)){
  for(k in 1:24){
    #i=658;k=5
    thispt <- summits[i,]
    
    geopt <- biomepts |> mutate(dst = ((thispt$lat-lat)^2+((thispt$lon-lon)*cos(thispt$lat*pi*2/360))^2)^0.5,
                                rnk = rank(dst), wts = case_when(rnk <= 10 ~ 100,
                                                                 rnk <= 20 ~ 1,
                                                                 TRUE ~ 0))|>st_drop_geometry()
    geopt$tval <- geopt[,thesecols[k]]
    mod <- lm(tval ~ lat+lon+elev, geopt, weights=geopt$wts)
    rate <- mod$coefficients[4]
   newelev <- thispt$summit
   oldelev <- thispt$elev
   oldtval <- st_drop_geometry((thispt[,thesecols[k]]))[,1]
    newtval <- (newelev-oldelev)*rate+oldtval
    summits2[i,thesecols[k]] <- newtval
  }}
summits2$elev <- summits2$summit

biomepts <- biomepts |> rbind(summits2[,colnames(biomepts)])
saveRDS(biomepts,'data/biomeptsummits.RDS')
#------------
library(climatools)
library(ggplot2)
library(sf)
library(terra)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

bpts <- readRDS('data/biomeptsummits.RDS')
bpts <- st_drop_geometry(bpts)
kth <- which(colnames(bpts) %in% 'th01')
ktl <- which(colnames(bpts) %in% 'tl01')
th <- bpts[,kth:(kth+11)]
tl <- bpts[,ktl:(ktl+11)]
kp <- which(colnames(bpts) %in% 'p01')
p <- bpts[,kp:(kp+11)]

tmean <- (th+tl)/2
bt <- sapply(tmean, function(x){ifelse(x > 0,x,0)})
tg <- pmax(apply(bt[,c(1:4,11:12)], MARGIN = 1, FUN = mean), apply(bt[,c(5:10)], MARGIN = 1, FUN = mean))
bpts$Tg <- tg
pet <- th
for(i in 1:12){#i=1
  pet[,i]<- GetPET(mon=i, lat=bpts$lat, th=th[,i], tl=tl[,i], p=p[,i])
}
petall <- apply(pet, MARGIN=1, FUN=sum)
pptall <- apply(p, MARGIN=1, FUN=sum)
Tcall <- apply(tmean, MARGIN=1, FUN=min)
Tclall <- apply(tl, MARGIN=1, FUN=min)
cindex <- pmin(Tcall,pmax(Tclall,Tcall-6)+3)

bpts <- readRDS('data/biomeptsummits.RDS')
bpts$Tg <- tg
bpts$m <- pptall/(pptall+petall)
bpts$cindex <- cindex
unique(st_drop_geometry(bpts[,c("BIOME","BIOME_NAME")]))







# ptssf3 <- bpts |> mutate(probiome = case_when(BIOME %in% c(1,3) & ((Tg >= 24 & elev < 1000) | (elev < 500 & Tg >= 12)) & 
#                                                   !(grepl('subtrop',tolower(ECO_NAME)) & lat > 22) ~ 'tropical moist forest',
#                                                 
#                                                 BIOME %in% c(1,4,5) & grepl('subtrop',tolower(ECO_NAME)) & lat >= 22 & elev < 500 ~ 'subtropical forest',
#                                                 BIOME %in% c(2) & Tg >= 6 ~ 'tropical dry forest',
#                                                 BIOME %in% c(1,3) & ((Tg < 24 & elev >= 500 & Tg >= 6)|(elev >= 1000 & Tg >= 6)) & !(grepl('subtrop',tolower(ECO_NAME)) & lat > 22) ~ 'tropical montane forest',
#                                                 BIOME %in% c(4,5) & elev < 500 & Tg >= 6 & grepl('coastal|valdiv|magell', tolower(ECO_NAME)) ~ 'coastal forest',
#                                                 BIOME %in% c(4,5) & (Tg >= 12 | (elev < 500 & Tg >= 6)) ~ 'temperate forest',
#                                                 
#                                                 BIOME %in% c(4,5) & Tg < 12 & Tg >= 6 & elev >= 500 ~ 'subalpine forest',
#                                                 BIOME %in% c(6) & Tg >= 6 ~ 'boreal forest',
#                                                 BIOME %in% c(11) & Tg < 6 & Tg > 0 ~ 'arctic tundra',
#                                                 BIOME %in% c(4,5,8,5) & Tg < 6 & elev >= 500 & Tg > 0 ~ 'alpine tundra',
#                                                 BIOME %in% c(1,2,3,7,10) & Tg < 6 & Tg > 0 & elev >= 1000 & lat < 22 & lat > -22 ~ 'paramo',
#                                                 BIOME %in% c(13) & m < 0.333  ~ 'desert',
#                                                 BIOME %in% c(7,8,10) & Tg >= 6  ~ 'grassland',
#                                                 BIOME %in% c(12) & Tg >= 6  ~ 'scrub',
#                                                 Tg == 0  ~ 'ice')) |> subset(!is.na(Tg) &  !is.na(Tg) & !is.na(probiome)) |> st_drop_geometry()

ptssf3 <- bpts |> mutate(veg1 = case_when(ECO_NAME %in% c('Guayanan Highlands moist forests','Florida sand pine scrub')|(grepl('scrub',ECO_NAME) & !grepl('forest',ECO_NAME)) ~ 'scrub',
                                          BIOME %in% c(1,2,3,4,5,6)  & m >= 0.333 ~ 'forest',
                                          BIOME %in% c(7,8,10) ~ 'grassland',
                                          BIOME %in% c(12,13) | (BIOME %in% c(1,2,3,4,5,6)  & m < 0.333) ~ 'scrub'),
                         winter1 = case_when(cindex < 0 & Tg-cindex >=12 ~ 'continental',
                                             cindex < 12 & Tg-cindex >=6 & Tg < 18 ~ 'oceanic',
                                             cindex < 12 & Tg-cindex >=6  ~ 'subtropical',
                                             Tg >= 18 | elev >= 1000 ~ 'tropical',
                                             TRUE ~ 'hyperoceanic'),
                         rain1=case_when(m < 0.333 ~ 'desert',
                                         m < 0.5 ~ 'dry',
                                         m < 0.667 | Tg < 6 | !veg1 %in% 'forest'~ 'moist',
                                         TRUE ~ 'rain'),
                         summer1=case_when(Tg < 6 & elev >= 1000 ~ 'alpine',
                                           Tg < 6  ~ 'arctic',
                                           Tg < 12  & elev >= 1000~ 'subalpine',
                                           Tg < 12   ~ 'boreal',
                                           Tg < 18  & elev >= 1000 ~ 'montane',
                                           Tg < 18   ~ 'hemiboreal',
                                           Tg < 24  & elev >= 1000 ~ 'submontane',
                                           Tg < 24   ~ 'warm',
                                           TRUE   ~ 'hot'),
                         probiome=case_when(Tg <= 0 ~ 'ice',
                                            Tg < 6 & winter1 %in% 'tropical' ~ trimws(paste(winter1,summer1,rain1,'paramo')),
                                            Tg < 6 ~ trimws(paste(winter1,summer1,rain1,'tundra')),
                                            m < 0.2 ~ trimws(paste(winter1,summer1,'desert')),
                                            !is.na(veg1) ~ trimws(paste(winter1,summer1,rain1,veg1))))
                         

                         unique(ptssf3$probiome)

                         
mon <- 1:12
mnn <- sprintf("%02d", mon) 
thesecols <- c(paste0('th',mnn),paste0('tl',mnn),paste0('p',mnn))
cnmsnum <- c('Tg','cindix','m')

od <- ptssf3 |> st_drop_geometry() |>  subset(!is.na(Tg) & !is.na(probiome)) |> group_by(winter1) |>
  mutate(cindex = median(cindex)) |> group_by(summer1) |> mutate(Tg = median(Tg)) |>
  group_by(rain1) |> mutate(m=median(m)) |> ungroup() |> subset(select=c(winter1, summer1, rain1, veg1, probiome, Tg, cindex, m)) |> unique()
od <- od |> arrange(cindex, m, Tg,veg1)
odprobiome <- unique(od$probiome)


cnmsnum <- c("lat","lon","elev","Tg","m","cindex",thesecols)
sumtab <- ptssf3 |> st_drop_geometry() |>  subset(!is.na(Tg) & !is.na(probiome)) |> 
  mutate(probiome=factor(probiome, levels = odprobiome)) |> group_by(LEVEL_NAME,REGION_NAM,CONTINENT,probiome) |>
  summarize(across(all_of(cnmsnum), median))|>ungroup() 

saveRDS(sumtab, 'data/biomesum.RDS')












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
      geom_point(data=sclim, mapping=aes(x=cindex, y=Tg), color = 'black', size=0.5)+
      #geom_density2d(data=StationMeans, mapping=aes(x=Cindex, y=Tg),color = 'black',alpha = 0.25)+
      scale_fill_manual("Legend", values = c("alpine" = "lightblue",
                                             "boreal" = "darkgreen",
                                             "continental" = "yellowgreen",
                                             "andean" = "pink",
                                             "oceanic" = "darkcyan",
                                             "subtropical" = "orange",
                                             "montane" = "pink4",
                                             "tropical" = "darkred"
      ))+
      
      scale_x_continuous(name= "Coldest Month (Daily Minimum)",
                         breaks=breaksx,
                         labels=paste0(breaksx,'(', breaksx-3,')'))+
      scale_y_continuous(name= "Growing Season", breaks=breaksy, expand=F)+
      coord_fixed(ratio = 1/1,xlim = c(-45,30), ylim = c(0, 33))+
      labs(title = paste("Climate of ",sname, "", sep=""))+
      theme_bw()+
      theme(legend.position='right',axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
            panel.grid.major = element_line(), panel.grid.minor = element_blank())
    
    
    #climplot3 ----
    #moisture x temperature
    bw1=data.frame(y=c(0,0,6,6), x=c(0,1,1,0))
    bw2=data.frame(y=c(6,6,12,12), x=c(0,1,1,0))
    bw3=data.frame(y=c(12,12,18,18), x=c(0,1,1,0))
    bw4=data.frame(y=c(18,18,24,24), x=c(0,1,1,0))
    bw5=data.frame(y=c(24,24,30,30), x=c(0,1,1,0))
    
    bmm1=data.frame(y=c(0,0,30,30), x=c(0,0.1111,0.1111,0))
    bmm2=data.frame(y=c(0,0,30,30), x=c(0.1111,0.2,0.2,0.1111))
    bmm3=data.frame(y=c(0,0,30,30), x=c(0.2,0.3333,0.333,0.2))
    bmm4=data.frame(y=c(0,0,30,30), x=c(0.3333,0.5,0.5,0.3333))
    bmm5=data.frame(y=c(0,0,30,30), x=c(0.5,0.6667,0.6667,0.5))
    bmm6=data.frame(y=c(0,0,30,30), x=c(0.6667,1,1,0.6667))
    
ggplot() +
      geom_polygon(data=bw1, mapping=aes(x=x, y=y, fill='alpine'),alpha = 0.2)+
      geom_polygon(data=bw2, mapping=aes(x=x, y=y, fill='cool'),alpha = 0.2)+
      geom_polygon(data=bw3, mapping=aes(x=x, y=y, fill='mild'),alpha = 0.2)+
      geom_polygon(data=bw4, mapping=aes(x=x, y=y, fill='warm'),alpha = 0.2)+
      geom_polygon(data=bw5, mapping=aes(x=x, y=y, fill='hot'),alpha = 0.2)+
      geom_polygon(data=bmm1, mapping=aes(x=x, y=y, fill='perarid'),alpha = 0.1)+
      geom_polygon(data=bmm2, mapping=aes(x=x, y=y, fill='arid'),alpha = 0.1)+
      geom_polygon(data=bmm3, mapping=aes(x=x, y=y, fill='semiarid'),alpha = 0.1)+
      geom_polygon(data=bmm4, mapping=aes(x=x, y=y, fill='subhumid'),alpha = 0.1)+
      geom_polygon(data=bmm5, mapping=aes(x=x, y=y, fill='humid'),alpha = 0.1)+
      geom_polygon(data=bmm6, mapping=aes(x=x, y=y, fill='perhumid'),alpha = 0.1)+
      geom_point(mapping=aes(x=m/100, y=Tg), color = 'black', shape="\u2605", size=3)+
      scale_fill_manual("Legend", values = c("alpine" = "cyan",
                                             "cool" = "green",
                                             "mild" = "yellow",
                                             "warm" = "orange",
                                             "hot" = "red",
                                             "perarid" = "red",
                                             "arid" = "orange",
                                             "semiarid" = "yellow",
                                             "subhumid" = "green",
                                             "humid" = "cyan",
                                             "perhumid" = "blue"
      ))+
      scale_y_continuous(name= "Growing Season", breaks=c(6,12,18,24,30),
                         labels=c('alpine/arctic 6', 'cool 12', 'mild 18','warm 24','hot 30'))+
      scale_x_continuous(name= "P/PET Ratio", breaks=c(0, .1111, .2,0.3333,0.5,0.6667),
                         labels=c('perarid', 'arid 0.125', 'semiarid 0.25','subhumid 0.5','humid 1','perhumid 2'))+
      coord_fixed(ratio = 1/30,ylim = c(0,30), xlim = c(0, 1))+
      
      labs(title = paste("Climate of ",sname, "", sep=""))+
      theme_bw()+
      theme(legend.position='none', axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1), 
            axis.text.y = element_text(vjust = 0), 
            panel.grid.major = element_line(), panel.grid.minor = element_blank()) 
    
    mtable <- df[,c('mon', 't', 'th', 'tl', 'p', 'e', 'a','u', 'Ra')]
    xtable::xtable(mtable)