# Loading required libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggimage)
library(magick)
library(grid)
library(shiny)
library(shinyWidgets)
library(tidyr)
library(readxl)

options(scipen = 80) # loading file size less than 50MB

file_pre <- "./" # same directory as shiny file


dt <- read_excel("shinydata1.xlsx", 
                 sheet = 1,      
                 col_names = TRUE,  
                 skip = 0,  
                 n_max = Inf)
                 # range = cell_cols(1:7))        

# 再进行pivot_longer
dt %>%
  pivot_longer(
    cols = -c(Organ, Age, Gender),
    names_to = "Metabolite", 
    values_to = "Intensity" 
  ) -> demo_df

Meta_v <- sort(unique(demo_df$Metabolite))

addResourcePath(
  prefix = "myimages",  
  directoryPath = paste0(file_pre,"www")
    )

ui <- fluidPage(
  
  # Application title
  titlePanel("The Organ-Resolved Metabolome Atlas of the Aging Mouse"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(width = 3,
                 h1("Parameter Setting", align = "center",
                    style = "font-family:arial; font-size:25pt; color:#000000"),
                 hr(),
                 pickerInput(
                   inputId = "meta",
                   label = "Metabolite",
                   choices = Meta_v,
                   multiple = FALSE,
                   options = list(
                     `actions-box` = TRUE,
                     size = 10),
                   selected = c("Carnitine")
                 ),
                 
                 p("Totally, ~1900 lipids, ~1000 metabolites.", 
                   style = "font-family: Arial; font-size: 14px;"),
                 
                 # br(),
                 hr(),
                 # self-adaption in image size 
                 p(img(src = "myimages/Overview.png", height = 200, width = "100%"), align = "left"),
                 # br(),
                 
                 p("A metabolome atlas for ~2900 structurally annotated metabolites 
                   of old versus young mice (90-96 weeks of age vs. 16 weeks) 
                   across 22 tissues and 4 biofluids.", 
                   style = "font-family: Arial; font-size: 14px;"),
                 
                 hr(),
                 awesomeCheckboxGroup(
                   inputId = "gender",
                   label = "Gender",
                   choices = c("Female", "Male"),
                   selected = c("Female", "Male"),
                   inline = TRUE
                 ),
                 hr()
    ),
    
    mainPanel(
      h1(textOutput("text1"), align = "center",
         style = "font-family:arial; font-size:30pt; color:#000000"), 
      hr(),
     
      
      fluidRow(
        column(width = 12, height = 550,
               plotOutput("plot_bar1", width = "auto", height = 450),
               br(), 
               hr()
        ),
        column(width = 12, height = 550,
               plotOutput("plot_bar2", width = "auto", height = 450),
               br(), 
               hr(), 
               # add comments
               p("1. Box-and-whisker plot, a statistical visualization tool, is used to summarize the distribution of metabolites intensity of organs.", 
                 style = "font-family:arial;font-size:12pt; color:#000000"),
               p("2. T-test was used for significance analysis. 
                     We use the following convention for symbols indicating statistical significance, ns: p > 0.05;  *: p <= 0.05;  **: p <= 0.01;  ***: p <= 0.001", 
                 style = "font-family:arial;font-size:12pt; color:#000000"),
               p("3. Explanation of Tissue Abbreviations: cerebrospinal fluid (CSF); brown adipose tissue (BAT); visceral adipose tissue (VAT); subcutaneous adipose tissue (SAT); cardiovascular (CV); Resp. (respiratory).", 
                 style = "font-family:arial;font-size:12pt; color:#000000")
        )
        
      )
    )
  )
)




fun1 <- function(df){

pathfig = file_pre

group1 <- c("Hippocampus",
            "CSF",
            "BAT",
            "SAT",
            "VAT",
            "Stomach",
            "Doudenum",
            "Jejunum",
            "Ileum",
            "Cecum",
            "Colon",
            "Liver",
            "Gallbladder"
)

group1_Fname <- c(
            "Hippocampus",
            "CSF",
            "BAT",
            "SAT",
            "VAT",
            "Stomach",
            "Duodenum",
            "Jejunum",
            "Ileum",
            "Cecum",
            "Colon",
            "Liver",
            "gallbladder"
)

group2 <- c("Spleen",
            "Thymus",
            "Testes",
            "Uterus",
            "Kidney",
            "Bladder",
            "Pancreas",
            "Heart",
            "Quadriceps",
            "Lung",
            "Plasma",
            "Urine",
            "Feces"
)

group2_Fname <- c(
            "Spleen",
            "Thymus",
            "Testes",
            "Uterus",
            "Kidney",
            "Bladder",
            "Pancreas",
            "Heart",
            "Quadriceps",
            "Lung",
            "Plasma",
            "Urine",
            "Feces"
)

image_grobs <- list()
for (i in seq_along(group1_Fname)) {
  img <- image_read((paste0(pathfig,"www/", group1_Fname[i],".png"))) %>%
    image_resize("600x") %>%
    image_rotate(0)
  grob <- rasterGrob(img, interpolate = TRUE)
  image_grobs[[i]] <- grob
}

image_grobs_2 <- list()
for (i in seq_along(group2_Fname)) {
  img <- image_read((paste0(pathfig,"www/", group2_Fname[i],".png"))) %>%
    image_resize("600x") %>%
    image_rotate(0)
  grob <- rasterGrob(img, interpolate = TRUE)
  image_grobs_2[[i]] <- grob
}

data2_group1 <- df %>%
  select(Organ, Age, Intensity) %>%
  filter(Organ %in% group1)

data2_group1$Organ <- factor(data2_group1$Organ,
                             levels = group1
)

data2_group1 %>%
  filter(Age == "Week 92") %>%
  group_by(Organ) %>%
  summarise(
    serine_mean = max(Intensity, na.rm = TRUE), 
    .groups = "drop"  
  ) %>%
  pull() %>%
  max() -> max1_group1

data2_group1 %>%
  filter(Age == "Week 16") %>%
  group_by(Organ) %>%
  summarise(
    serine_mean = max(Intensity, na.rm = TRUE), 
    .groups = "drop"  
  ) %>%
  pull() %>%
  max() -> max2_group1

p_values_group1 <- data2_group1 %>%
  group_by(Organ) %>%
  do({
    group_week16 <- .$Intensity[.$Age == "Week 16"]
    group_week92 <- .$Intensity[.$Age == "Week 92"]
    
      if(var(group_week16, na.rm = TRUE) == 0 & var(group_week92, na.rm = TRUE) == 0) {
      data.frame(p = NA_real_)
    } else {
      t_test_result <- t.test(Intensity ~ Age, data = .)
      data.frame(p = t_test_result$p.value)
    }
  }) %>%
  mutate(
    p_label = case_when(
      is.na(p) ~ "ns",  
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    y_pos = max(max1_group1, max2_group1)  
  )

anno_data_group1 <- data.frame(
  x = c(1.5, 4, 9), 
  y = rep(p_values_group1$y_pos[1] , 3), 
  label = c("Nervous", "Connective Tissue", "Digestive")
)


data2_group2 <- df %>%
  select(Organ, Age, Intensity) %>%
  filter(Organ %in% group2)

data2_group2$Organ <- factor(data2_group2$Organ,
                             levels = group2
)

data2_group2 %>%
  filter(Age == "Week 92") %>%
  group_by(Organ) %>%
  summarise(
    serine_mean = max(Intensity, na.rm = TRUE), 
    .groups = "drop"  
  ) %>%
  pull() %>%
  max() -> max1_group2

data2_group2 %>%
  filter(Age == "Week 16") %>%
  group_by(Organ) %>%
  summarise(
    serine_mean = max(Intensity, na.rm = TRUE), 
    .groups = "drop"  
  ) %>%
  pull() %>%
  max() -> max2_group2


p_values_group2 <- data2_group2 %>%
  group_by(Organ) %>%
  # 先检查每组数据是否有变异
  do({
    # 提取两组数据
    group_week16 <- .$Intensity[.$Age == "Week 16"]
    group_week92 <- .$Intensity[.$Age == "Week 92"]
    
    # 检查两组数据是否都为恒定值
    if(var(group_week16, na.rm = TRUE) == 0 & var(group_week92, na.rm = TRUE) == 0) {
      # 两组都无变异，返回p=NA
      data.frame(p = NA_real_)
    } else {
      # 至少一组有变异，进行t检验
      t_test_result <- t.test(Intensity ~ Age, data = .)
      data.frame(p = t_test_result$p.value)
    }
  }) %>%
  mutate(
    p_label = case_when(
      is.na(p) ~ "ns",  # 无变异时标记为ns
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    y_pos = max(max1_group2, max2_group2)  
  )

anno_data_group2 <- data.frame(
  x = c(1.5, 3.5, 5.5, 7, 8, 9, 10, 12), 
  y = rep(p_values_group2$y_pos[1], 8), 
  label = c("Immune", "Reproductive", "Excretory",
            "Endocrine",
            "CV",
            "Muscular",
            "Resp.",
            "Biofluid & digestive product")
)

return(list(data2_group1, p_values_group1, anno_data_group1, image_grobs,
            data2_group2, p_values_group2, anno_data_group2, image_grobs_2))

}


server <- function(input, output) {
  reaction1 <- reactive({
    req(input$meta, input$gender) # the trick
    options(scipen = 80)
    demo_df %>%
      filter(Metabolite == input$meta,
             Gender %in% input$gender) %>%
      fun1()
  })

  output$text1 <- renderText({
    input$meta
  })
  
    output$plot_bar1 <- renderPlot({
    
    ggboxplot(reaction1()[[1]],
              x = "Organ",
              y = "Intensity",
              bxp.errorbar = T,
              fill = "Age",
              palette = "npg",
              border.size = 1,
              add.params = list(
                size = 2, 
                width = 0.3 
              )
    ) +
      geom_text(
        data = reaction1()[[2]],
        aes(x = Organ, y = y_pos*1.05, label = p_label),
        position = position_dodge(0.8), 
        size = 5
      ) +
        scale_y_continuous(
          breaks = function(x) pretty(x, n = 6),
          labels = scales::label_scientific()
        ) +
      theme(
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(angle = 45, size = 12, hjust = 1), 
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)
      ) +
      labs(x = "Organ",
           y = "Intensity") +
      geom_segment(
        data = reaction1()[[2]],
        aes(x = 1, xend = 2,
            y = reaction1()[[2]]$y_pos[1] *1.3,
            yend = reaction1()[[2]]$y_pos[1]*1.3),
        size = 0.2
      ) +
      geom_segment(
        data = reaction1()[[2]],
        aes(x = 3, xend = 5,
            y = reaction1()[[2]]$y_pos[1] *1.3,
            yend = reaction1()[[2]]$y_pos[1]*1.3),
        size = 0.2
      ) +
      geom_segment(
        data = reaction1()[[2]],
        aes(x = 6, xend = 13,
            y = reaction1()[[2]]$y_pos[1] *1.3,
            yend = reaction1()[[2]]$y_pos[1]*1.3),
        size = 0.2
      ) +
      geom_text(
        data = reaction1()[[3]],
        aes(x = x, y = y*1.4, label = label),
        size = 5
      ) +
      annotation_custom(
        grob = reaction1()[[4]][[1]],
        xmin = 0.8,
        xmax = 1.2,
        ymin = reaction1()[[2]]$y_pos[1] *1,
        ymax = reaction1()[[2]]$y_pos[1] *1.3
      ) +
      annotation_custom(
        grob = reaction1()[[4]][[2]],
        xmin = 1.8,
        xmax = 2.2,
        ymin = reaction1()[[2]]$y_pos[1] *1,
        ymax = reaction1()[[2]]$y_pos[1] *1.3
      ) +
      annotation_custom(
        grob = reaction1()[[4]][[3]],
        xmin = 2.8,
        xmax = 3.2,
        ymin = reaction1()[[2]]$y_pos[1] *1,
        ymax = reaction1()[[2]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[4]][[4]],
        xmin = 3.8,
        xmax = 4.2,
        ymin = reaction1()[[2]]$y_pos[1] *1,
        ymax = reaction1()[[2]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[4]][[5]],
        xmin = 4.8,
        xmax = 5.2,
        ymin = reaction1()[[2]]$y_pos[1] *1,
        ymax = reaction1()[[2]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[4]][[6]],
        xmin = 5.8,
        xmax = 6.2,
        ymin = reaction1()[[2]]$y_pos[1] *1,
        ymax = reaction1()[[2]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[4]][[7]],
        xmin = 6.8,
        xmax = 7.2,
        ymin = reaction1()[[2]]$y_pos[1] *1,
        ymax = reaction1()[[2]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[4]][[8]],
        xmin = 7.8,
        xmax = 8.2,
        ymin = reaction1()[[2]]$y_pos[1] *1,
        ymax = reaction1()[[2]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[4]][[9]],
        xmin = 8.8,
        xmax = 9.2,
        ymin = reaction1()[[2]]$y_pos[1] *1,
        ymax = reaction1()[[2]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[4]][[10]],
        xmin = 9.8,
        xmax = 10.2,
        ymin = reaction1()[[2]]$y_pos[1] *1,
        ymax = reaction1()[[2]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[4]][[11]],
        xmin = 10.8,
        xmax = 11.2,
        ymin = reaction1()[[2]]$y_pos[1] *1,
        ymax = reaction1()[[2]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[4]][[12]],
        xmin = 11.8,
        xmax = 12.2,
        ymin = reaction1()[[2]]$y_pos[1] *1,
        ymax = reaction1()[[2]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[4]][[13]],
        xmin = 12.8,
        xmax = 13.2,
        ymin = reaction1()[[2]]$y_pos[1] *1,
        ymax = reaction1()[[2]]$y_pos[1] *1.3
      )

  })
  
  
  output$plot_bar2 <- renderPlot({
    
    ggboxplot(reaction1()[[5]],
              x = "Organ",
              y = "Intensity",
              bxp.errorbar = T,
              fill = "Age",
              palette = "npg",
              border.size = 1,
              add.params = list(
                size = 2, 
                width = 0.3 
              )
    ) +
      geom_text(
        data = reaction1()[[6]],
        aes(x = Organ, y = y_pos*1.05, label = p_label),
        position = position_dodge(0.8), 
        size = 5
      ) +
      scale_y_continuous(
        breaks = function(x) pretty(x, n = 6),
        labels = scales::label_scientific()
      ) +
      theme(
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(angle = 45, size = 12, hjust = 1),  
        axis.text.y = element_text(size = 12),  
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt")
      ) +
      labs(x = "Organ",
           y = "Intensity") +
      geom_segment(
        data = reaction1()[[6]],
        aes(x = 1, xend = 2, 
            y = reaction1()[[6]]$y_pos[1] *1.3, 
            yend = reaction1()[[6]]$y_pos[1]*1.3),
        size = 0.2  
      ) +
      geom_segment(
        data = reaction1()[[6]],
        aes(x = 3, xend = 4, 
            y = reaction1()[[6]]$y_pos[1] *1.3, 
            yend = reaction1()[[6]]$y_pos[1]*1.3),
        size = 0.2 
      ) +
      geom_segment(
        data = reaction1()[[6]],
        aes(x = 5, xend = 6, 
            y = reaction1()[[6]]$y_pos[1] *1.3, 
            yend = reaction1()[[6]]$y_pos[1]*1.3),
        size = 0.2 
      ) +
      geom_segment(
        data = reaction1()[[6]],
        aes(x = 6.8, xend = 7.2, 
            y = reaction1()[[6]]$y_pos[1] *1.3, 
            yend = reaction1()[[6]]$y_pos[1]*1.3),
        size = 0.2 
      )  +
      geom_segment(
        data = reaction1()[[6]],
        aes(x = 7.8, xend = 8.2, 
            y = reaction1()[[6]]$y_pos[1] *1.3, 
            yend = reaction1()[[6]]$y_pos[1]*1.3),
        size = 0.2
      )  +
      geom_segment(
        data = reaction1()[[6]],
        aes(x = 8.8, xend = 9.2, 
            y = reaction1()[[6]]$y_pos[1] *1.3, 
            yend = reaction1()[[6]]$y_pos[1]*1.3),
        size = 0.2 
      )  +
      geom_segment(
        data = reaction1()[[6]],
        aes(x = 9.8, xend =10.2, 
            y = reaction1()[[6]]$y_pos[1] *1.3, 
            yend = reaction1()[[6]]$y_pos[1]*1.3),
        size = 0.2 
      )  +
      geom_segment(
        data = reaction1()[[6]],
        aes(x = 10.8, xend = 13.2, 
            y = reaction1()[[6]]$y_pos[1] *1.3, 
            yend = reaction1()[[6]]$y_pos[1]*1.3),
        size = 0.2 
      ) +
      geom_text(
        data = reaction1()[[7]],
        aes(x = x, y = y*1.4, label = label),
        size = 5
      ) +
      annotation_custom(
        grob = reaction1()[[8]][[1]],
        xmin = 0.8,
        xmax = 1.2,
        ymin = reaction1()[[6]]$y_pos[1] *1,
        ymax = reaction1()[[6]]$y_pos[1] *1.3
      ) +
      annotation_custom(
        grob = reaction1()[[8]][[2]],
        xmin = 1.8,
        xmax = 2.2,
        ymin = reaction1()[[6]]$y_pos[1] *1,
        ymax = reaction1()[[6]]$y_pos[1] *1.3
      ) +
      annotation_custom(
        grob = reaction1()[[8]][[3]],
        xmin = 2.8,
        xmax = 3.2,
        ymin = reaction1()[[6]]$y_pos[1] *1,
        ymax = reaction1()[[6]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[8]][[4]],
        xmin = 3.8,
        xmax = 4.2,
        ymin = reaction1()[[6]]$y_pos[1] *1,
        ymax = reaction1()[[6]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[8]][[5]],
        xmin = 4.8,
        xmax = 5.2,
        ymin = reaction1()[[6]]$y_pos[1] *1,
        ymax = reaction1()[[6]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[8]][[6]],
        xmin = 5.8,
        xmax = 6.2,
        ymin = reaction1()[[6]]$y_pos[1] *1,
        ymax = reaction1()[[6]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[8]][[7]],
        xmin = 6.8,
        xmax = 7.2,
        ymin = reaction1()[[6]]$y_pos[1] *1,
        ymax = reaction1()[[6]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[8]][[8]],
        xmin = 7.8,
        xmax = 8.2,
        ymin = reaction1()[[6]]$y_pos[1] *1,
        ymax = reaction1()[[6]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[8]][[9]],
        xmin = 8.8,
        xmax = 9.2,
        ymin = reaction1()[[6]]$y_pos[1] *1,
        ymax = reaction1()[[6]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[8]][[10]],
        xmin = 9.8,
        xmax = 10.2,
        ymin = reaction1()[[6]]$y_pos[1] *1,
        ymax = reaction1()[[6]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[8]][[11]],
        xmin = 10.8,
        xmax = 11.2,
        ymin = reaction1()[[6]]$y_pos[1] *1,
        ymax = reaction1()[[6]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[8]][[12]],
        xmin = 11.8,
        xmax = 12.2,
        ymin = reaction1()[[6]]$y_pos[1] *1,
        ymax = reaction1()[[6]]$y_pos[1] *1.3
      )+
      annotation_custom(
        grob = reaction1()[[8]][[13]],
        xmin = 12.8,
        xmax = 13.2,
        ymin = reaction1()[[6]]$y_pos[1] *1,
        ymax = reaction1()[[6]]$y_pos[1] *1.3
      ) 
    
  })
  
}


shinyApp(ui = ui, server = server)
