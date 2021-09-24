shinyAppServer <- function(input, output) {

  # create reactive value to store variables
  react_val = reactiveValues()

  observe({

    if (!is.null(input$config_file)) {

      source(input$config_file$datapath)

      #check input directory
      if (!dir.exists(input_dir)) {
        input_dir <- "not_found"
      }

      #check output directory
      if (!dir.exists(output_dir)) {
        output_dir <- "not_found"
      }

      #check reference file
      if (!file.exists(reference)) {
        reference <- "not_found"
      }

      #check annotation directory
      if (!file.exists(annotation_dir)) {
        annotation_dir <- "not_found"
      }

      #check reference version
      if (!reference_version %in% c("hg19", "hg38")) {
        reference_version <- "not_found"
      }

      react_val$input_dir <- input_dir
      react_val$output_dir <- output_dir
      react_val$reference <- reference
      react_val$annotation_dir <- annotation_dir
      react_val$reference_version <- reference_version

    }
  })

  output$config_table <- renderTable(colnames = FALSE, rownames = TRUE, {

    if(!is.null(input$config_file)) {
      config_table <- data.frame(values = c(react_val$input_dir, react_val$output_dir, react_val$reference, react_val$annotation_dir))
      rownames(config_table) <- c("Input directory", "Output directory", "Genome reference", "SNP annotation directory")
      config_table
    }

  })

  output$memory <- renderText({
    paste0("Mutect2 will require approximately: ", input$parallel*15, "GB of RAM \nTotal system RAM is approximately: ", round(as.numeric(grep("[0-9]+", strsplit(system("cat /proc/meminfo", intern = TRUE)[1], split = " ")[[1]], value = TRUE))/1000000))
  })


  output$metadata_table <- renderTable({

    if (is.null(input$metadata_file)) {
      return(NULL)
    } else {
      react_val$metadata = read.table(input$metadata_file$datapath, header = TRUE)
      head(react_val$metadata)
    }

  })

  observeEvent(input$SNPsom_wrapper,  {
    SNPsom_wrapper(react_val$metadata, react_val$input_dir, react_val$output_dir, react_val$reference, input$threads, input$parallel, react_val$reference_version, react_val$annotation_dir)
  })

}
