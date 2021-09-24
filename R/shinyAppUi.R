shinyAppUi <- fluidPage(

  titlePanel("SNPsom"),

  sidebarLayout(

    sidebarPanel(

      fileInput("config_file", label = h4("Config file")),
      br(),
      fileInput("metadata_file", label = h4("Metadata file")),
      br(),
      numericInput("threads", label = h4("Number of threads"), min = 1, max = grep("[0-9]+", strsplit(system("lscpu | egrep '^CPU\\(s\\)'", intern = TRUE), split = " ")[[1]], value = TRUE), value = 2),
      br(),
      numericInput("parallel", label = h4("Number of parallel samples for Mutect2"), min = 1, max = floor((as.numeric(grep("[0-9]+", strsplit(system("cat /proc/meminfo", intern = TRUE)[1], split = " ")[[1]], value = TRUE))/1000000)/15), value = 2),
      textOutput(outputId = "memory"),
      br(),
      actionButton("SNPsom_wrapper", label = "Start SNPsom analysis")

    ),

    mainPanel(

      tableOutput(outputId =  "config_table"),

      tableOutput(outputId = "metadata_table")

    )

  )
)
