# AKSTA Case Study 1 (Group 2)
# authors Tess Landon (01617344), Simona Hovančíková(12347335), Babak Bayani (12347302)

library(shiny)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)
library(countrycode)
library(jsonlite)
library(tidyr)

# Load data
data_cia<-fromJSON("data_cia2.json")
data_cia<-as.data.frame(data_cia)

# Rename columns
data_cia<-data_cia %>% rename(
    Country = country,Continent = continent,
    `Education Expenditure (% GDP)`=expenditure,
    `Youth Unemployment Rate (%)`=youth_unempl_rate,
    `Net Migration Rate (%)`= net_migr_rate,
    `Electricity from Fossil Fuels (%)`=electricity_fossil_fuel,
    `Population Growth Rate (%)`=pop_growth_rate,
    `Life Expectancy`=life_expectancy,
    Population=population,Area=area)

variables<-c(
  "Education Expenditure (% GDP)",
  "Youth Unemployment Rate (%)",
  "Net Migration Rate (%)",
  "Electricity from Fossil Fuels (%)",
  "Population Growth Rate (%)",
  "Life Expectancy")

ui<-fluidPage(
  titlePanel(div(
    h1("CIA World Factbook 2020 Insights"), h5(em('Welcome! This Shiny App - created by yours Truly - allows you to explore information from the CIA 2020 Factbook through various visualizations, descriptive uni- or multivariate statistics, and other features. Enjoy! ')))),
  tabsetPanel(tabPanel("Univariate Analysis",
                sidebarLayout(
                 sidebarPanel(
                   selectInput("uni_var", "Select a variable:", choices = variables),
                   actionButton("view_data", "View raw data"),
                   DTOutput("raw_data")),
                 mainPanel(tabsetPanel(
                   tabPanel("Map", plotlyOutput("map_plot")),
                   tabPanel("Global Analysis", plotlyOutput("global_boxplot"), plotlyOutput("global_density")),
                   tabPanel("Analysis per Continent", plotlyOutput("continent_boxplot"), plotlyOutput("continent_density"))
                   )))),
      tabPanel("Multivariate Analysis",
               sidebarLayout(
                 sidebarPanel(
                   selectInput("var_x", "Select X variable:", choices = variables),
                   selectInput("var_y", "Select Y variable:", choices = variables),
                   selectInput("size_var", "Size by:", choices = c("Population", "Area"))),
                 mainPanel(plotlyOutput("scatter_plot"))))))

server <- function(input,output){
  selected_data <- reactive({
    req(input$uni_var)
    data_cia %>% select(Country,Continent,all_of(input$uni_var))})
  
  observeEvent(input$view_data,{
    output$raw_data<-renderDT({
      datatable(head(selected_data(), nrow(data_cia)))})})
  
  output$map_plot <- renderPlotly({
    world_map <- map_data("world")
    world_map$ISO3 <- countrycode(world_map$region, "country.name", "iso3c")
    
    map_data_joined <- left_join(world_map, data_cia, by = c("ISO3" = "ISO3"))
    p <- ggplot(map_data_joined, aes(x = long, y = lat, group = group, fill = .data[[input$uni_var]])) +
      geom_polygon(color = "white") +
      scale_fill_viridis_c(option = "C", na.value = "grey50") +theme_void()
    ggplotly(p, tooltip = c("region", input$uni_var))})
  
  output$global_boxplot <- renderPlotly({
    p <- ggplot(data_cia, aes(y =.data[[input$uni_var]])) +
      geom_boxplot(fill = "lightblue") +
      theme_minimal()
    ggplotly(p)})
  
  output$global_density<-renderPlotly({
    p<-ggplot(data_cia, aes(x = .data[[input$uni_var]]))+
      geom_histogram(aes(y=..density..),bins=30,fill="skyblue",alpha=0.6)+
      geom_density(alpha = 0.5, fill = "darkblue")+theme_minimal()
    ggplotly(p)})
  
  output$continent_boxplot<-renderPlotly({
    p<-ggplot(data_cia,aes(x=Continent,y=.data[[input$uni_var]],fill=Continent))+
      geom_boxplot()+theme_minimal()
    ggplotly(p)})
  
  output$continent_density<-renderPlotly({
    p <- ggplot(data_cia, aes(x =.data[[input$uni_var]], fill = Continent)) +
      geom_density(alpha = 0.5) +
      theme_minimal()
    ggplotly(p)})
  
  output$scatter_plot<-renderPlotly({
    p<-ggplot(data_cia,aes(x=.data[[input$var_x]],y=.data[[input$var_y]],color=Continent,
                           size=.data[[input$size_var]],group=Continent))+geom_point(alpha=0.6)+
      geom_smooth(method="loess",se=FALSE,aes(color=Continent,group=Continent))+theme_minimal()
    ggplotly(p,tooltip=c("label", "x", "y", "color"))})}

# Run the Shiny app
shinyApp(ui = ui, server = server)
