#Install required libraries (uncomment if you need to install any package)
# install.packages(c("shiny", "ggplot2", "ggpubr", "rstan", "dplyr",
#                    "tidyr","reshape2", "writexl" ))

# Load required libraries
library(shiny)
library(ggplot2)
library(ggpubr)
library(rstan)
library(dplyr)
library(tidyr)
library(reshape2)
library("writexl")

# Define UI
ui <- fluidPage(
  titlePanel("Dengue FOI and Model Analysis"),
  sidebarLayout(
    sidebarPanel(
      # Input: File selector for the first CSV file
      fileInput("dengue_cases", "Choose the Dengue Cases CSV file",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      
      # Input: File selector for the second CSV file
      fileInput("population", "Choose the Population CSV file",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      
      # Radio buttons for selecting the separator used in the CSV file
      # radioButtons("sep", "Separator",
      #              choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
      #              selected = ","),
      
      # Dynamic input for selecting locations (Districts)
      # uiOutput("location_select"),
      uiOutput("province_select"),
      
      # Dynamic input for selecting years
      uiOutput("year_select"),
      
      # Dynamic input for reference age groups
      uiOutput("reference_agegroup_select"),
      
      # Model selection
      selectInput("model", "Select Model", choices = c("time_varying.stan", "time_constant.stan")),
      
      sliderInput("n_iter",
                  "Number of iterations to run:",
                  min = 5000,
                  max = 10000,
                  value = 5000,
                  step = 1000),
      
      # Button to run the model
      actionButton("runModel", "Run Model")
      # downloadButton("downloadAllPlots", "Download All Plots")
    ),
    
    mainPanel(
      uiOutput("rhatWarning"),  # <- Barra di avviso condizionale
      tabsetPanel(
        # tabPanel("Model Feedback", verbatimTextOutput("modelStatus")),  # New feedback tab
        tabPanel("Fit Plot", plotOutput("fitPlot")),
        tabPanel("FOI Plot", plotOutput("lambdaPlot")),
        tabPanel("Pars Plot", plotOutput("parsPlot")),
        tabPanel("Model Feedback", verbatimTextOutput("modelStatus")),  # New feedback tab
        tabPanel("Traceplot", plotOutput("tracePlot")),  # New tab for Traceplot
        tabPanel("Rhat Values", tableOutput("rhatTable"))  # New tab for Rhat values
      )
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  
  # Helper function to read CSV file
  read_file <- function(file) {
    if (is.null(file)) return(NULL)
    read.csv(file$datapath, sep = ",", stringsAsFactors = FALSE)
  }
  
  # Reactive expressions for data
  dengue_cases_data <- reactive({
    read_file(input$dengue_cases)
  })
  
  population_data <- reactive({
    read_file(input$population)
  })
  
  # Dynamically update the location input based on uploaded Dengue Cases file
  # output$location_select <- renderUI({
  #   req(dengue_cases_data())
  #   districts <- unique(dengue_cases_data()$District)
  #   selectInput("location", "Select Location(s)", 
  #               choices = districts, 
  #               selected = districts[1], 
  #               multiple = TRUE)
  # })
  
  output$province_select <- renderUI({
    req(dengue_cases_data())
    
    # Get unique provinces
    provinces <- unique(dengue_cases_data()$Province)
    
    # Create a list where each province is displayed along with its districts
    province_choices <- lapply(provinces, function(province) {
      # Get districts belonging to the current province
      districts <- dengue_cases_data() %>%
        filter(Province == province) %>%
        pull(District) %>%
        unique()
      
      # Create a label for each province with its districts listed
      paste(province, "(", paste(districts, collapse = ", "), ")")
    })
    
    # Use selectInput to display provinces with associated districts in labels
    selectInput("province", "Select Province", choices = setNames(provinces, province_choices), selected = provinces[1], multiple = TRUE)
  })
  
  
  # Dynamically update the year input based on uploaded Dengue Cases file
  output$year_select <- renderUI({
    req(dengue_cases_data())
    
    years <- sort(unique(dengue_cases_data()$Year))  # order years
    
    # sliderfor time window selection
    sliderInput("year", "Select Year Range", 
                min = min(years), 
                max = max(years), 
                value = c(min(years), max(years)),  # initial value is across all years
                step = 1,                           # each step is 1 years
                sep = "")                           # remove comma 
  })
  
  output$reference_agegroup_select <- renderUI({
    req(dengue_cases_data())
    
    agegroups <- unique(dengue_cases_data()$AgeGroup)  # order years
    agegroups <- agegroups[-c(1, length(agegroups))] #extreme excluded
    # sliderfor time window selection
    
    selectInput("Ref_Agegroup", "Select Reference Age-group", 
                choices = agegroups, 
                selected = agegroups[4], 
                multiple = F)
  })
  
  
  # Reactive expressions for filtered data
  filtered_cases <- reactive({
    req(dengue_cases_data())
    dengue_cases_data() %>%
      filter(Province %in% input$province, Year %in% c(seq(input$year[1],input$year[2],1)))
  })
  
  filtered_population <- reactive({
    req(population_data())
    population_data() %>%
      filter(Province %in% input$province, Year %in% c(seq(input$year[1],input$year[2],1)))
  })
  
  # Run the model when the button is clicked
  model_fit <- eventReactive(input$runModel, {
    req(filtered_cases())
    req(filtered_population())
    
    
  
    
    # plot_ready(FALSE)  # Hide plots
    
    # Model Feedback Output: Reset the status output
    output$modelStatus <- renderText({"Model is starting..."})
    
    withProgress(message = 'Fitting Model', value = 0, {
      # withProgress(message = 'Fitting Model', value = 0, {
      n_steps <- 5
      #   
      # Increment progress for data preparation
      incProgress(1/n_steps, detail = "Preparing data...")
      
      # Data preparation steps (as shown in your original code)
      df_long_cases <- filtered_cases() %>%
        select(-c("X","X.1"))
      df_long_pop <- filtered_population()  %>%
        select(-c("X","X.1"))
      
      age_group <- unique(df_long_cases$AgeGroup)
      age_min <-  unique(df_long_cases$min_age) +1 
      age_max <- unique(df_long_cases$max_age) +1
      
      df_long <- merge(df_long_cases, df_long_pop, by = c("District", "Year", "Province", "AgeGroup"))
      df_long$Cases <- as.numeric(df_long$Cases)
      df_long$Population <- as.numeric(df_long$Population)
      
      province_list <- unique(df_long$Province)
      
      df_long <- df_long %>%
        filter(!is.na(Cases), !is.na(Population), Population > 0)
      
      df_long$incidence <- as.integer((df_long$Cases / df_long$Population) * 10000)
      
      output$modelStatus <- renderText({"Incidence calculated"})
      
      # Prepare cases matrix for the model
      case_wide <- df_long %>%
        select( -"min_age.x", -"max_age.x", -"min_age.y", - "max_age.y" ,  -"Province", - "Population", -"Cases") %>%  # Exclude these columns
        pivot_wider(names_from = "AgeGroup", values_from = "incidence")
      
      df_long <- melt(case_wide, id.vars = c("District", "Year"))
      cases <- acast(df_long, District ~ Year ~ variable, value.var = "value")
      
      #associate province number to each district (dinamic dictionary based on selected provinces)
      selected_province_dict <-setNames(seq_along(unique(df_long_cases$Province)), unique(df_long_cases$Province))
      
      # Extract district order from cases dataframe
      districts_order <- dimnames(cases)[[1]] 
      
      output$modelStatus <- renderText({"prob2"})
      
      # Create province vector assocaited to districts
      province_numbers <- sapply(districts_order, function(d) {
        province <- df_long_cases$Province[df_long_cases$District == d][1]
        if (province %in% input$province) {
          selected_province_dict[[province]]  
        } else {
          NA  # it should not happen
        }
      })
      
      
      # Increment progress after data preparation
      # incProgress(1/n_steps, detail = "Data prepared, running Stan model...")
      
      nP <- length(unique(input$province))
      nA <- dim(cases)[3]
      nT <- length(unique(df_long$Year))
      nL <- length(unique(df_long$District))
      
      age_band_width <- age_max - age_min + 1
      ref <- grep(input$Ref_Agegroup, age_group)
      
      stan_data <- list(
        nA = nA, nL = nL, nT = nT, 
        age = seq(0, 99, 1), 
        X = as.numeric(ref),
        # N = as.numeric(input$N),
        nP = nP,
        N = 10000,
        prov_N = province_numbers,
        cases_inc = cases, 
        ageLims = rbind(age_min, age_max),
        sum_cases = apply(cases, c(1, 2), sum),
        age_band = age_band_width
      )
      
      
      output$modelStatus <- renderText({"Data list ready, model is starting..."})
      
      # Running the Stan model
      fit <- sampling(
        stan_model(file = input$model),
        data = stan_data,
        seed = 12345,  # Set seed here
        chains = 3, iter = input$n_iter, warmup = 1000, cores = 3, refresh = 100
      )
      
      # Increment progress after Stan model run
      incProgress(2/n_steps, detail = "Processing Stan output...")
      
      # Process the results (as in your original code)
      chains <- rstan::extract(fit)
      
      
      #extrat historical lambda estimates
      lam_H <- as.data.frame(matrix(NA, nL, 3))
      chi1 <- as.data.frame(matrix(NA, nL, 3))
      chi2 <- as.data.frame(matrix(NA, nL, 3))
      rho <- as.data.frame(matrix(NA, nP, 3))
      
      
      for(l in 1:nL) {
        lam_H[l,1:3] <- quantile(chains$lam_H[,l],c(0.5,0.025,0.975))
        lam_H <- as.data.frame(lam_H)
        colnames(lam_H)[1:3] <- c("med", "ciL", "ciU")
        lam_H$loc[l] <- districts_order[l]
        lam_H$year <- "lam_H"
        
        chi1[l,1:3] <- quantile(chains$chi1[,l],c(0.5,0.025,0.975))
        chi1 <- as.data.frame(chi1)
        colnames(chi1)[1:3] <- c("med", "ciL", "ciU")
        chi1$loc[l] <- districts_order[l]
        chi1$year <- "chi1"
        
        chi2[l,1:3] <- quantile(chains$chi2[,l],c(0.5,0.025,0.975))
        chi2 <- as.data.frame(chi2)
        colnames(chi2)[1:3] <- c("med", "ciL", "ciU")
        chi2$loc[l] <- districts_order[l]
        chi2$year <- "chi2"
        
      }
      
      for(l in 1:nP){
        rho[l,1:3] <- quantile(chains$rho[,l],c(0.5,0.025,0.975))
        rho <- as.data.frame(rho)
        colnames(rho)[1:3] <- c("med", "ciL", "ciU")
        rho$loc[l] <- province_list[l]
        rho$year <- "rho"
      }
      
      
      gamma <- as.data.frame(matrix(NA, nP, 3))
      for (l in 1:nP) {
        gamma[l, 1:3] <- quantile(chains$gamma[, l], c(0.5, 0.025, 0.975))
        gamma <- as.data.frame(gamma)
        colnames(gamma)[1:3] <- c("med", "ciL", "ciU")
        gamma$loc[l] <- province_list[l]
        gamma$year <- "gamma"
      }
      
      pars <- rbind(rho, gamma, lam_H, chi1, chi2)
      
      # Select model-specific parameters
      if (input$model == "time_varying.stan") {
        lam <- array(NA, dim = c(nT, 5, nL))
        
        for (l in 1:nL) for (t in 1:nT) {
          lam[t, 4, l] <- districts_order[l]
          lam[t, 5, l] <- unique(df_long$Year)[t]
          lam[t, 1:3, l] <- quantile(chains$lam_t[, l, t], c(0.5, 0.025, 0.975))
        }
        
        lambda <- do.call(rbind, lapply(1:nL, function(l) lam[,, l]))
        lambda <- as.data.frame(lambda)
        colnames(lambda) <- c("med", "ciL", "ciU", "loc", "year")
        lambda$year <- rep(seq(input$year[1], input$year[2], 1), nL)
        
        lambda$med <- as.numeric(lambda$med)
        lambda$ciL <- as.numeric(lambda$ciL)
        lambda$ciU <- as.numeric(lambda$ciU)
        
        writexl::write_xlsx(lambda, path = "lambda_estimates.xlsx")
      } else {
        
        lambda <- lam_H
      }
      
      # Write lambda estimates to Excel
      
      
      writexl::write_xlsx(pars, path = "SL_pars_estimates.xlsx")
      
      output$modelStatus <- renderText({
        paste("Model fitting completed successfully at", Sys.time())
      })
      
      # plot_ready(TRUE)  # Print plots
      
      return(list(fit = fit, lambda = lambda, case_wide = cases, pars= pars, 
                  nT=nT, nL=nL, nA=nA, district_list = districts_order, age_group = age_group))
    })
  })
  
  
  # define tables for lambda and pars
  
  output$lambdaTable <- renderTable({
    if (input$model == "time_varying.stan"){
      lambda_values <- as.data.frame(model_fit()$lambda)
    } else  {
      lambda_values <- as.data.frame(model_fit()$lambda) 
    }
  })
  
  output$parsTable <- renderTable({
    pars_values <- as.data.frame(model_fit()$pars)
  })
  
  # Plot Outputs
  output$fitPlot <- renderPlot({
    
    fit <- model_fit()
    chains <- rstan::extract(fit$fit)
    
    fit_tot <- data.frame()
    
    if (input$model == "time_varying.stan"){
      
      for (i in 1:length(model_fit()$district_list)) {
        fit.df <- data.frame(model_fit()$case_wide[i,,])
        fit.df <- fit.df %>%
          pivot_longer(c(1:model_fit()$nA), names_to = "Age_Group", values_to = "Cases")
        
        fit.df$year <- rep(1:model_fit()$nT, each = model_fit()$nA)
        
        ageG <- fit.df$Age_Group[1:model_fit()$nA]
        fit.df[, c('pred', 'ciL', 'ciU')] <- NA
        
        fit.df <- as.data.frame(fit.df)
        for (t in 1:model_fit()$nT) for (a in 1:model_fit()$nA) {
          fit.df[fit.df$year == t & fit.df$Age_Group == ageG[a], 4:6] <- quantile(chains$Ecases[, i, t, a], c(0.5, 0.025, 0.975))
        }
        
        fit.df$Age_group <- rep(model_fit()$age_group, model_fit()$nT)
        # fit.df$Age_group <- factor(fit.df$Age_group, levels = c("0-04 yrs", "05-09 yrs", "10-14 yrs", "15-19 yrs", "20-24 yrs", "25-49 yrs", "50-64 yrs", "65+ yrs"))
        
        fit.df$loc <- unique(model_fit()$district_list)[i]
        fit_tot <- rbind(fit_tot, fit.df)
      }
      
      fit_tot$year <- fit_tot$year 
      
      
      
      plot <-  ggplot(fit_tot) +
        geom_point(aes(x = Age_group, y = Cases, colour = "observed"), size = 2) +
        geom_ribbon(aes(x = Age_group, ymin=ciL, ymax=ciU, group=year, fill= "estimated 95% CrI"), alpha = 0.3) +
        geom_line(aes(x = Age_group, y = pred, group=year, colour = "estimated median"), size = 1) +
        facet_grid(loc ~ year, scales = "free_y") +
        labs(x = NULL, y = paste("Incidence per 10.000 population")) +
        scale_color_manual(name = NULL, values = c("observed" = "black", "estimated median" = "dodgerblue")) +
        scale_fill_manual(name = NULL, values = c("observed" = "black", "estimated 95% CrI" = "dodgerblue")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top"
        )
      
    } else {
      
      for (i in 1:length(model_fit()$district_list)) {
        fit.df <- data.frame(model_fit()$case_wide[i,,])
        fit.df <- fit.df %>%
          pivot_longer(c(1:model_fit()$nA), names_to = "Age_Group", values_to = "Cases")
        
        fit.df$year <- rep(1:model_fit()$nT, each = model_fit()$nA)
        
        ageG <- fit.df$Age_Group[1:model_fit()$nA]
        fit.df <- as.data.frame(fit.df)
        
        fit.df <- fit.df %>%
          group_by(Age_Group) %>%
          summarise(Cases = mean(Cases), .groups ="drop")
        
        fit.df$year <- 1
        fit.df[, c('pred', 'ciL', 'ciU')] <- NA
        
        for (a in 1:model_fit()$nA) {
          fit.df[fit.df$Age_Group == ageG[a], "pred"] <- quantile(chains$Ecases[, i, a], 0.5)
          fit.df[fit.df$Age_Group == ageG[a], "ciL"] <- quantile(chains$Ecases[, i, a], 0.025)
          fit.df[fit.df$Age_Group == ageG[a], "ciU"] <- quantile(chains$Ecases[, i, a], 0.975)
          # fit.df[fit.df$Age_Group == ageG[a], c("pred", "ciL", "ciU")] <- list(quantile(chains$Ecases[, i, a], c(0.5, 0.025, 0.975)))
          
        }
        
        fit.df$Age_group <- c("0-04 yrs", "05-09 yrs", "10-14 yrs", "15-19 yrs", "20-24 yrs", "25-49 yrs", "50-64 yrs", "65+ yrs")
        fit.df$Age_group <- factor(fit.df$Age_group, levels = c("0-04 yrs", "05-09 yrs", "10-14 yrs", "15-19 yrs", "20-24 yrs", "25-49 yrs", "50-64 yrs", "65+ yrs"))
        
        fit.df$loc <- unique(model_fit()$district_list)[i]
        fit_tot <- rbind(fit_tot, fit.df)
      }
      
      # fit_tot$year <- fit_tot$year + 2016
      plot <- ggplot(fit_tot)+
        geom_point(aes(x = Age_group, y = Cases, colour = "observed"), size = 2) +
        geom_line(aes( x = Age_group, y = pred, group=year, colour= "estimated median"), size = 1)+
        geom_ribbon(aes(x = Age_group, ymin=ciL, ymax=ciU, group=year, fill= "estimated 95% CrI"),alpha=0.3)+
        facet_wrap(~loc, ncol =  model_fit()$nT, scales = "free_y") +
        labs(x = NULL, y = paste("Incidence per 10.000 population")) +
        scale_color_manual(name = NULL, values = c("observed" = "black", "estimated median" = "dodgerblue")) +
        scale_fill_manual(name = NULL, values = c("observed" = "black", "estimated 95% CrI" = "dodgerblue")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top"
        )
      
      
    }
    plot
    
    return(plot)
  }, height = function() {
    800 + 20 * length(unique(model_fit()$nL))  # dynamic height
  })
  
  output$lambdaPlot <- renderPlot({
    
    
    if (input$model == "time_varying.stan"){
      lambda <- model_fit()$lambda
      lambda.years <- unique(lambda$year)
      lamH <- filter(model_fit()$pars, year == "lam_H")
      lambda <- rbind(lambda, lamH)
      lambda$year <- factor(lambda$year, levels = c("lam_H", lambda.years))
    } else  {
      lambda <- model_fit()$lambda
    }
    
    
    plot <- ggplot(lambda, aes(x = year, y = med, col = year)) +
      geom_point() +
      geom_errorbar(aes(ymin = ciL, ymax = ciU), width = 0.2) +
      labs(title = NULL, x = NULL, y = NULL, col = NULL) +
      # scale_color_manual(values = colors_scale) +
      facet_wrap(~ loc)+
      theme_bw() +
      # ylim(0,0.1)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none")
    plot
    
    return(plot)
  })
  
  output$parsPlot <- renderPlot({
    pars <- model_fit()$pars
    pars1 <- filter(pars, year %in% c("rho", "gamma"))
    pars2 <- filter(pars, year %in% c("chi1", "chi2"))
    
    plot1 <- ggplot(pars1, aes(x = loc, y = med, col = year)) +
      geom_point() +
      geom_errorbar(aes(ymin = ciL, ymax = ciU), width = 0.2) +
      labs(title = NULL, x = NULL, y = NULL, col = NULL) +
      # scale_color_manual(values = colors_scale) +
      facet_wrap(~ year, ncol = 2, scale = "free_y")+
      theme_bw() +
      ylim(0,1)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none")
    
    plot2 <- ggplot(pars2, aes(x = loc, y = med, col = year)) +
      geom_point() +
      geom_errorbar(aes(ymin = ciL, ymax = ciU), width = 0.2) +
      labs(title = NULL, x = NULL, y = NULL, col = NULL) +
      # scale_color_manual(values = colors_scale) +
      facet_wrap(~ year, ncol=2, scale = "free_y")+
      theme_bw() +
      # ylim(0,1)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none")
    
    plot <- ggarrange(plot1, plot2, ncol=1)
    return(plot)
  })
  
  # Traceplot Output
  output$tracePlot <- renderPlot({
    fit <- model_fit()$fit
    
    # Extract the traceplot for specific parameters
    traceplot(fit, inc_warmup = FALSE, ncol = 4)
    
  })
  
  # Rhat Values Output
  output$rhatTable <- renderTable({
    fit <- model_fit()$fit
    
    # Extract Rhat values from the model
    fit_summary <- summary(fit)$summary
    rhat_values <- fit_summary[, "Rhat", drop = FALSE]
    
    # Convert to a data frame for better display
    rhat_df <- as.data.frame(rhat_values)
    rhat_df <- tibble::rownames_to_column(rhat_df, var = "Parameter")
    
    rhat_df  # Return Rhat values as a table
  })
  
  output$rhatWarning <- renderUI({
    req(model_fit())  # Assicurati che il modello sia stato eseguito
    
    # Estrai Rhat dal fit
    fit <- model_fit()$fit
    fit_summary <- summary(fit)$summary
    max_rhat <- max(fit_summary[, "Rhat"], na.rm = TRUE)
    
    if (max_rhat > 1.2) {
      # Barra di avviso rossa
      div(
        style = "background-color: #FFCCCC; color: red; padding: 10px; margin-bottom: 10px; border: 1px solid red; border-radius: 5px;",
        strong("WARNING: Improve model fit (Rhat > 1.2)! suggestion: change reference age group and increase iterations")
      )
    } else {
      NULL  # Nessun messaggio se Rhat va bene
    }
  })
  
  # Server
  output$downloadAllPlots <- downloadHandler(
    filename = function() {
      paste("Dengue_Model_Plots_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      # Open a PDF device to save plots
      pdf(file, width = 11, height = 8.5)  # Landscape layout
      
      # Generate each plot and print it to the PDF
      print(output$fitPlot())
      print(output$lambdaPlot())
      print(output$parsPlot())
      print(output$tracePlot())
      
      # Close the PDF device after saving all plots
      dev.off()
    }
  )
  
}


# Run the application
shinyApp(ui = ui, server = server)


print("DONE")
