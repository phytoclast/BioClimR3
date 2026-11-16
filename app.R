library(shiny)
library(climatools)
library(ggplot2)
biosum <- readRDS('data/biomesum.RDS')
ui <- bslib::page_fluid(
  
  tags$head(
    tags$style(HTML("

      body { font-size: 12px; }
    "))
  ),
  
  headerPanel('Climate Workshop'),
  sidebarLayout(
    sidebarPanel(
      fluidRow(  
        column(10,htmlOutput("country"))),
      fluidRow(    
        column(10,htmlOutput("biome"))),
      fluidRow(
        radioButtons(
        'soil','soil texture', choices = c('gravel'=0.03,sand=0.06535,'loam'=0.1104667,'muck'=0.5413333), selected =0.06535),
        sliderInput(
        'root', 'Effective root depth (cm)', 100, min = 0, max = 500, step = 10)),
      fluidRow(
        verbatimTextOutput('AWC'))
    ),
    mainPanel(
      fluidRow(
        column(12,plotOutput("climplot",  width = "100%", height = "600px"))),
      fluidRow(
        column(12,verbatimTextOutput('climatext')))
      )))



server <- function(input, output) {
  rv <- reactiveValues(my_text = "") 
  output$country = renderUI({
    
    listcountries <- sort(unique(biosum$LEVEL_NAME))
    
    selectInput(inputId = "countryselect", #name of input
                label = "Select Name:", #label displayed in ui
                choices = unique(listcountries), #calls list of available counties
                selected = 'Michigan')
  })
  output$biome = renderUI({
    
    listbiomes1 <- subset(biosum, LEVEL_NAME %in% input$countryselect)
    listbiomes <- sort(unique(listbiomes1$probiome))
    selectInput(inputId = "biomeselect", #name of input
                label = "Select Biome:", #label displayed in ui
                choices = unique(listbiomes), #calls list of available counties
                selected = 'temperate forest')
  })
  
  output$AWC <- renderText({
    paste("Available water holding capacity",round(as.numeric(input$soil)*input$root*10,1),"mm.")
  })
  
  output$climplot <- renderPlot({
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
    
    awc <- as.numeric(input$soil)*input$root*10
    
    df$Ra <- GetSolarRad(GetDayNumber(df$mon),Lat=lat)
    df$e <- GetPET.df(mon=df$mon,lat=lat,th=df$th,tl=df$tl,p=df$p)
    df<-df|>cbind(SoilWaterBudget(df$p, df$e, awc=awc))
    df$a <- pmin(df$e, df$p)
    Tg <- climatools::warm6ppt(t=df$t,p=0,n=6,max_var = TRUE, v_input = 't', v_return = 'b',bt=0)
    MAAT <- mean(df$t)
    Tp <- weighted.mean(df$t, df$p)
    Tw <- max(df$t)
    Tc <- min(df$t)
    Tcl <- min(df$tl)
    Twh <- max(df$th)
    PET <- sum(df$e)
    MAP <- sum(df$p)
    AET <- sum(df$a)
    prat <- MAP/(PET+0.01)
    m <- MAP/(PET+MAP+0.01)*100
    d <- 100-AET/(PET+0.01)*100
    Tclx <- climatools::XtremLow(Tcl=Tcl,Lat=lat,Lon=lon,Elev=elev)
    cindex <- pmin(Tc,Tcl+3)
    cindex2 <- pmin(Tc, Tclx*12/15+12)
    trindex <- Tg-cindex
    climatext <- paste0("lat: ",round(lat,1),"° | lon: ",round(lon,1),"° | elev: ",round(elev,0)," m| years: 1961-1990\n",
                        "MAAT: ",round(MAAT,1),"°C | MAP: ",round(MAP,0)," mm | PET: ",round(PET,0)," mm\n",
                        "Twh: ",round(Twh,1),"°C | Tw: ",round(Tw,1),"°C | Tg: ",round(Tg,1),"°C\n",
                        "Tc: ",round(Tc,1),"°C | Tcl: ",round(Tcl,1),"°C | Tclx: ",round(Tclx,1),"°C\n",
                        "Tp: ",round(Tp,1),"°C | MAP/PET: ",round(prat,2)," | moisture index: ",round(m,0),"% | drought index: ",round(d,0),"%")
    
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
    
    
    rv$climatext <- climatext
    
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
  })
  
  output$climatext <- renderText({
    rv$climatext
    
  })
}
shinyApp(ui, server)
