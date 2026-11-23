library(shiny)
library(climatools)
library(ggplot2)
library(sf)
biosum <- readRDS('data/biomesum.RDS')
lakes <- readRDS(file='data/lakes.RDS')
states <- readRDS(file='data/states.RDS')

ui <- bslib::page_fluid(
  
  tags$head(
    tags$style(HTML("

      body { font-size: 12px; }
    "))
  ),
  
  headerPanel('Climate Workshop'),
  sidebarLayout(
    sidebarPanel(
      ###############################
      fluidRow(  
        column(10,htmlOutput("country"))),
      fluidRow(    
        column(10,htmlOutput("biome"))),
      ###############################
      fluidRow(
        radioButtons(
          'soil','soil texture', choices = c('gravel'=0.04,sand=0.06535,'loam'=0.1104667,'muck'=0.5413333), selected =0.06535),
        sliderInput(
          'root', 'Effective root depth (cm)', 100, min = 0, max = 500, step = 10)),
      fluidRow(
        verbatimTextOutput('AWC')),
      fluidRow(
        column(12,verbatimTextOutput('climatext'))),
      fluidRow(
        radioButtons("RadioGraphtype",inline = T,  label = ("Select Graph"),
                     choiceNames = list(HTML("<font size=-3>Monthly"), 
                                        HTML("Summer × Winter"), 
                                        HTML("Summer × Moisture"), 
                                        HTML("Map")
                     ),
                     choiceValues = list(1,2,3,4),
                     selected = 1)),
      HTML("</font>")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel('Plot',
                 fluidRow(column(12,plotOutput("climplot",  width = "100%", height = "600px"))),
                 fluidRow(column(width = 2,checkboxInput("simplebar", label = "Simple Bar Graph", value = FALSE)))),

        tabPanel('Table',tableOutput('Climtab')),
        tabPanel('Reference',verbatimTextOutput('climaref'))
      )
    )))



server <- function(input, output) {
  rv <- reactiveValues(my_text = "") 
  ###############################
  output$country = renderUI({
    
    listcountries <- sort(unique(biosum$LEVEL_NAME))
    
    selectInput(inputId = "countryselect", #name of input
                label = "Select Name:", #label displayed in ui
                choices = unique(listcountries), #calls list of available countries
                selected = 'Alaska')
  })
  output$biome = renderUI({
    
    listbiomes1 <- subset(biosum, LEVEL_NAME %in% input$countryselect)
    listbiomes <- sort(unique(listbiomes1$probiome))
    selectInput(inputId = "biomeselect", #name of input
                label = "Select Biome:", #label displayed in ui
                choices = unique(listbiomes), #calls list of available countries
                selected = 'continental boreal dry forest')
  })
  ###############################
  output$AWC <- renderText({
    paste("Available water holding capacity",round(as.numeric(input$soil)*input$root*10,1),"mm.")
  })
  
  output$climplot <- renderPlot({
    ###############################
    sclim <- subset(biosum, LEVEL_NAME %in% input$countryselect & probiome %in% input$biomeselect)[1,]
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
    ###############################
    awc <- as.numeric(input$soil)*input$root*10
    
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
    pAET <- sum(df$aa)
    AET <- sum(df$a)
    prat <- MAP/(PET+0.01)
    m <- MAP/(PET+MAP+0.01)*100
    d <- 100-pAET/(PET+0.01)*100
    Sd <- 100-AET/(PET+0.01)*100
    # Tclx <- climatools::XtremLow(Tcl=Tcl,Tc=Tc,Tw=Tw)
    Tclx <- climatools::XtremLow(Tcl=Tcl,Lat=lat,Lon=lon,Elev=elev)
    # cindex <- pmin(Tc,Tcl+3)
    cindex <- pmin(Tc,pmax(Tcl,Tc-6)+3)
    # cindex <- pmin(Tc, Tclx*12/15+12)
    trindex <- Tg-cindex
    climatext <- paste0("lat: ",round(lat,1),"° | lon: ",round(lon,1),"° | elev: ",round(elev,0)," m| years: 1961-1990\n",
                        "MAAT: ",round(MAAT,1),"°C | MAP: ",round(MAP,0)," mm | PET: ",round(PET,0)," mm\n",
                        "Twh: ",round(Twh,1),"°C | Tw: ",round(Tw,1),"°C | Tg: ",round(Tg,1),"°C\n",
                        "Tc: ",round(Tc,1),"°C | Tcl: ",round(Tcl,1),"°C | Tclx: ",round(Tclx,1),"°C\n",
                        "Tp: ",round(Tp,1),"°C | MAP/PET: ",round(prat,2)," | m: ",round(m,0),"% | d: ",round(d,0),"% | Sd: ",round(Sd,0)," mm")
    
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
    
    tindex <- case_when(cregime %in% 1 ~ paste(tindex,'Continental'),
                        cregime %in% 2 & Tg < 18 ~ paste(tindex,'Oceanic'),
                        cregime %in% 2 ~ paste(tindex,'Subtropical'),
                        cregime %in% 3 ~ paste(tindex,'Tropical'))
    
    
    alt <- case_when(cregime  %in% c(1,2) & Tg < 6 ~ 'Arctic/Alpine',
                     cregime  %in% 1 & Tg < 12 ~ 'Boreal/Subalpine',
                     cregime  %in% 1  ~ 'Continental',
                     cregime  %in% 2 & Tg < 18 ~ 'Oceanic',
                     cregime  %in% 2 ~ 'Subtropical',
                     cregime  %in% 3 & Tg < 6 ~ 'Andean',
                     cregime  %in% 3 & Tg < 24  ~ 'Tropical Montane',
                     cregime  %in% 3  ~ 'Tropical')
    
    #koppen
    df$rnk <- rank(-df$t, ties.method='first')
    nmonth10 <- max(subset(df, t >= 10)$rnk)
    pptsummer <- sum(subset(df,(mon %in% c(4,5,6,7,8,9) & lat >=0) | (mon %in% c(10,11,12,1,2,3) & lat < 0))$p)
    psummer <- pptsummer/(MAP+0.01) 
    pthresh <- MAAT*20+case_when(psummer >= 0.7 ~ 280,
                                 psummer >= 0.3 ~ 140,
                                 TRUE ~ 0)
    psmax <- max(subset(df,(mon %in% c(4,5,6,7,8,9) & lat >=0) | (mon %in% c(10,11,12,1,2,3) & lat < 0))$p)
    psmin <- min(subset(df,(mon %in% c(4,5,6,7,8,9) & lat >=0) | (mon %in% c(10,11,12,1,2,3) & lat < 0))$p)
    pwmax <- max(subset(df,(mon %in% c(10,11,12,1,2,3) & lat >=0) | (mon %in% c(4,5,6,7,8,9) & lat < 0))$p)
    pwmin <- min(subset(df,(mon %in% c(10,11,12,1,2,3) & lat >=0) | (mon %in% c(4,5,6,7,8,9) & lat < 0))$p)
    
    k1 <- case_when(nmonth10 < 1 ~ 'E',
                    MAP < pthresh ~ 'B',
                    Tc < 0 ~ 'D',
                    Tc < 18 ~ 'C',
                    TRUE  ~ 'A')
    
    
    k2 <- case_when(k1 %in% 'E' & Tw < 0 ~ 'F',
                    k1 %in% 'E'  ~ 'T',
                    k1 %in% 'B' & MAP < pthresh/2  ~ 'W',
                    k1 %in% 'B'  ~ 'S',
                    k1 %in% c('C','D') & psmin < 30 & psmin < pwmax/3  ~ 's',
                    k1 %in% c('C','D') & pwmin < psmax/10   ~ 'w',
                    k1 %in% c('C','D') ~ 'f',
                    pmin(pwmin,psmin) >= 60 ~ 'f',
                    pmin(pwmin,psmin) >= 100-MAP/25 ~ 'm',
                    TRUE ~ 'w')
    k3 <- case_when(k1 %in% c('E','A') ~ '',
                    k1 %in% 'B' & MAAT < 18   ~ 'k',
                    k1 %in% 'B'  ~ 'h',
                    Tc <= -38 ~ 'd',
                    k1 %in% c('C','D') & nmonth10 < 4 ~ 'c',
                    k1 %in% c('C','D') & Tw < 22 ~ 'b',
                    TRUE ~ 'a')
    
    climaclass <- paste(sname,'\n',mindex, tindex, paste0(' ',k1,k2,k3))
    climaref <- paste(
      "#Summary Statistics#\n
      MAAT = mean annual air temperature.\n
    MAP = mean annual precipitation.\n
    PET = potential evapotranspiration (annual).\n
    Twh = mean daily high temperature for the warmest month.\n
    Tw = mean temperature for the warmest month.\n
    Tc = mean temperature for the coldest month.\n
    Tg = growing season temperature (mean positive temperature for the warmest 6 months).\n
    Tcl = mean daily low temperature for the coldest month.\n
    Tclx = mean annual extreme low temperature.\n
  Tp = precipitation weighted mean temperature.\n
   m = moisture index, MAP/(PET+MAP). \n
   d = drought index, the total monthly deficit in precipitation relative to potential evapotranspiration.\n
   Sd = seasonal soil moisture deficit.\n\n
   #Monthly Table#\n
  mon = month.\n
  t = monthly mean temperature.\n
  th = monthly mean daily high temperature.\n
  tl = monthly mean daily low temperature.\n
  p = monthly mean precipitation.\n
  e = monthly PET.\n
  a = monthly actual evapotranspiration (considers both precipitation and soil water).\n
  u = monthly precipitation surplus.\n
  Ra = mean daily solar radiation.")
    
    rv$climaref <-climaref
    
    rv$climatext <- climatext
    
    climplot1 <- makeClimplot( t=df$t,
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
                               bar = input$simplebar)
    
    a1=data.frame(x=c(-50,-50,0,-6), y=c(0,6,6,0))
    a2=data.frame(x=c(-50,-50,0,-6), y=c(6,12,12,6))
    a3=data.frame(x=c(-50,-50,0,0), y=c(12,36,36,12))
    a4=data.frame(x=c(-6,0,6,0), y=c(0,6,6,0))
    a5=data.frame(x=c(-6,0,0,12,0), y=c(6,12,18,18,6))
    a6=data.frame(x=c(0,0,12,12), y=c(18,36,36,18))
    a7=data.frame(x=c(0,12,12,24,18,6), y=c(6,18,24,24,18,6))
    a8=data.frame(x=c(12,12,36,24), y=c(24,36,36,24))
    breaksx=seq(-60,30,6)
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
    
    climplot2 <-  
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
      geom_point(mapping=aes(x=cindex, y=Tg), color = 'black', size=5, shape =16)+
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
      scale_y_continuous(name= "Growing Season", breaks=breaksy, expand=FALSE)+
      coord_fixed(ratio = 1/1,xlim = c(-45,30), ylim = c(0, 33))+
      labs(title = paste("Climate of ",sname, "", sep=""))+
      theme_bw()+
      theme(legend.position='right',axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
            panel.grid.major = element_line(), panel.grid.minor = element_blank())
    
    
    
    #######-----
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
    
  climplot3 <- ggplot() +
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
    geom_point(mapping=aes(x=m/100, y=Tg), color = 'black', size=5, shape =16)+
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
    scale_y_reverse(name= "Growing Season", breaks=c(0,6,12,18,24,30),
                       labels=c('ice 0','alpine/arctic 6', 'cool 12', 'mild 18','warm 24','hot 30'), expand=FALSE)+
    scale_x_continuous(name= "P/PET Ratio", breaks=c(0, .1111, .2,0.3333,0.5,0.6667),
                       labels=c('perarid', 'arid 0.125', 'semiarid 0.25','subhumid 0.5','humid 1','perhumid 2'), expand=FALSE)+
    coord_fixed(ratio = 1/30,ylim = c(0,30), xlim = c(0, 1))+
    
    labs(title = paste("Climate of ",sname, "", sep=""))+
    theme_bw()+
    theme(legend.position='none', axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1), 
          axis.text.y = element_text(vjust = 0), 
          panel.grid.major = element_line(), panel.grid.minor = element_blank()) 
    #######-----
  graphymax = (sclim$lat)+10# 
  graphymin = (sclim$lat)-10# 
  graphxmax = (sclim$lon)+10/cos(sclim$lat*2*pi/360)# 
  graphxmin = (sclim$lon)-10/cos(sclim$lat*2*pi/360)# 
  ocean = data.frame(y=c(-90,-90,90,90), x=c(-180,180,180,-180))
  
  climmap <-  ggplot() +
    geom_polygon(data = ocean,
                 aes(x = x, y = y),
                 color = 'darkgray', fill = 'lightcyan', size = .2)+
    geom_sf(data = states,
            color = 'darkgray', fill = 'lightyellow', size = .2)+
    geom_sf(data = lakes,
            color = 'darkgray', fill = 'lightcyan', size = .2)+
    geom_point(data=sclim, mapping=aes(x=lon, y=lat), color = 'red', size=2)+

    coord_sf(xlim = c(graphxmin,graphxmax), ylim = c(graphymin,graphymax)) + theme_void()
  #######-----
  
    mtable <- df[,c('mon', 't', 'th', 'tl', 'p', 'e', 'a','u', 'Ra')]
    rv$my_table <- xtable::xtable(mtable)
    if(input$RadioGraphtype == 1){
            climplot1
    }else if(input$RadioGraphtype == 2) {
      climplot2
    }else if(input$RadioGraphtype == 3) {
      climplot3
    }else if(input$RadioGraphtype == 4) {
      climmap
    }
    
    
    
    
  })
  
  output$climatext <- renderText({
    rv$climatext
    
  })
  output$climaref <- renderText({
    rv$climaref
    
  })
  output$Climtab = renderTable({ 
    rv$my_table
  })
}
shinyApp(ui, server)
