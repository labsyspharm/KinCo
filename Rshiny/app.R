library(shiny)
library(r3dmol)
library(bio3d)
library(colourpicker)
library(stringr)
library(shinyWidgets)
library(shinythemes)
library(shinyBS)
library(stringr)
library(filesstrings)
library(R.utils)
library(reticulate)
library(dplyr)
library(stringr)

#WHEN PUBLISHED ON SHINY IO
py_install(c('pandas')) 
#LOCALLY
#use_python("/opt/anaconda3/envs/Py38.torch/bin/python")

source('molViewer.R')
source('export_posefile.R')

pd <- import('pandas')
#gdown <- import("gdown")
zipfile <- import("zipfile")

#UI elements
mainUI<-function(id){ #This is main UI
  
  
  
  ns <- NS(id)
  
    
    navbarPage(theme = shinytheme("cosmo"),
               title = "KinCo",
               
               selected = "3D Viewer", collapsible = F, inverse = T,
               tags$head(tags$style(".navbar {background-color:black;}")),
               tabPanel(div( tags$script(HTML("var header = $('.navbar > .container-fluid');
                header.append('<div style=\"float:right\"><img src=\"HMS.png\" alt=\"alt\" style=\"float:right;width:40px;height:50px;padding-top:0px;\"> </a>`</div>');
                console.log(header)")
               ),
               tags$script(HTML("var header = $('.navbar > .container-fluid');
                header.append('<div style=\"float:right\"><img src=\"CU.png\" alt=\"alt\" style=\"float:right;width:40px;height:50px;padding-top:0px;\"> </a>`</div>');
                console.log(header)")
               ),
               tags$script(HTML("var header = $('.navbar > .container-fluid');
                header.append('<div style=\"float:right\"><img src=\"Logo_3.png\" alt=\"alt\" style=\"float:right;width:40px;height:50px;padding-top:0px;\"> </a>`</div>');
                console.log(header)")
               ))),
               #tags$div(img(src='HMS.png',style="margin-top: -14px; padding-right:10px;padding-bottom:10px", height = 60)),
               tabPanel("About"),
               tabPanel("3D Viewer",
                        # sidebarPanel(
                        #   fileInput("file", "File input:"),
                        #   textInput("txt", "Text input:", "general"),
                        #   sliderInput("slider", "Slider input:", 1, 100, 30),
                        #   tags$h5("Default actionButton:"),
                        #   actionButton("action", "Search"),
                        #   
                        #   tags$h5("actionButton with CSS class:"),
                        #   actionButton("action2", "Action button", class = "btn-primary")
                        # ),
                        tags$h3(strong("Welcome to KinCo!")),
                        tags$h5("A database for", tags$em('in silico'), "kinase-compound structures"),
                        #tags$h5("Instructions to use KinCo:"),
                        br(),
                        br(),
                        fluidPage(
                          fluidRow(column(3,
                                          selectizeInput(  inputId=ns("kinaseSelector"),
                                                           #label = paste0("Choose a kinase to start with ",icon("caret-circle-down")),
                                                           label = span(tagList("Choose a kinase to start with  "," ",icon("arrow-right"))),
                                                           choices= kinases,
                                                           selected = NULL,
                                                           multiple = FALSE,
                                                        
                                                           options = list(
                                                             onInitialize = I('function() { this.setValue(""); }')
                                                          #   plugins = list("remove_button")
                                                           )
                                          )),
                                   column(1),
                                   column(4,
                                          selectizeInput(inputId=ns("compoundSelector"),
                                                         #label = paste0("Choose a kinase to start with ",icon("caret-circle-down")),
                                                         label = span(tagList("Next, choose a compound  "," ",icon("arrow-right"))),
                                                         choices= NULL,
                                                         selected = NULL,
                                                         multiple = FALSE,)
                          ),
                          column(1),
                          column(3,
                                 div(
                                   style="width:230px;height:30%",fluidRow(
                                     tags$h5(strong("Kinase-compound pair affinity (nM)")),
                                   verbatimTextOutput(ns("affinity"), placeholder = F)))
                        
                                 
                                         
                                    
                                  )
                          
                                   ),
                          br(),
                           fluidRow(column(2,
                                           actionButton(ns("submit"), "Generate Homologs?"))
                             
                           ),
                          br(),
                          br(),
                          fluidRow(column(12,
                                          div(id = ns('kinasePlaceholder'))
                          ))
                        )
               ),
               tabPanel("Funding"),
               tabPanel("Download")
               )
  
}

####server code kinase UI
kinaseUi <- function(id){
  ns <- NS(id)
  fluidPage(
    fluidRow(tags$h5(strong('Generating homolog models for the following pair: '))),
    
    fluidRow(tags$h5(textOutput(ns("pairChosen")))),
    #fluidRow(tags$h5('Locating our database... ')),
    fluidRow(tags$h5(textOutput(ns("download")))),
    #fluidRow(tags$h5('Please be patient...')),
    br(),
    br(),
    fluidRow(column(3,
                    
                    actionButton(inputId=ns("addItem"), "Add New Homolog")),
            #column(),
    fluidRow(column(3,
                    actionButton(inputId=ns("removeItem"), "Remove Homolog")))),
    br(),
    br(),
    fluidRow(column(4,
                    div(id = ns('homologModulePlaceholder_first'))),
    column(4,
           div(id = ns('homologModulePlaceholder_second'))),
    column(4,
           div(id = ns('homologModulePlaceholder_third')))
    )
  )
  
}

kinaseServer<-  function(id, kinase_compound_pair){
  
  
  
  moduleServer(
    id,
    function(input, output, session) {
      
      #posedf <- NULL
      
      counter<-reactiveValues()
      
      counter$count=0
      
      ns <-session$ns
      
      
      
      output$pairChosen <- renderText({paste0(kinase_compound_pair[1], ' - ', kinase_compound_pair[2])})
      
      compound_file = paste0(kinase_compound_pair[2],'.zip')
      compoundMetaData = kinco_meta[(kinco_meta$INCHIKEY_DESALT == compound_file) & (kinco_meta$Symbol == kinase_compound_pair[1]),]
      #homolog_id = kinco_meta[(kinco_meta$INCHIKEY_DESALT == compound_file)& (kinco_meta$Symbol == kinase_compound_pair[1]),]$PDB_GID
      
      
      showModal(modalDialog("Locating our database.......", footer=NULL))
      ########Retrieving information on cloud###########
      
      GENEID <<- compoundMetaData$ENTREZ_GENE_ID
      tgtgeneid <<- compoundMetaData$ENTREZ_GENE_ID
      lgdinchikey <<- kinase_compound_pair[2]
      datasetdir <<- './data/' 
      #pairpath <<- paste0(datasetdir, tgtgeneid, '/', lgdinchikey)
      pairpath <<- paste0(datasetdir,lgdinchikey,"/")
      output_dir <<- paste0(datasetdir,'exported_structures/')
      
      if (!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      ###pairdir <<- unzip_lgd(pairpath)
      
      
      
      
      unlink(paste0("./data/",GENEID))
      dir.create(paste0("./data/",GENEID))
      ###gdown$download(paste0("https://drive.google.com/uc?id=",compoundMetaData$GID,"&export=download"), output = paste0("./data/",GENEID,"/",compound_file))
      unzip("KinCo_sel_2022-08-14.zip", files = c(paste0("KinCo_sel_2022-08-14/",GENEID,".zip")), exdir = "data")
      unzip(paste0("data/KinCo_sel_2022-08-14/",GENEID,".zip"), files = c(compound_file), exdir = paste0("./data/",GENEID))
      unzip(paste0("data/",GENEID,"/",compound_file), exdir = paste0("./data/",lgdinchikey))
      posedf <<- pd$read_pickle(paste0(pairpath,'pose.df'))
      unlink( list.files("data", pattern = "\\.zip$", full.names = T))
      
      showModal(modalDialog(HTML("Locating our database... <br>
                            Kinase found! <br>
                            Searching for compound..."), footer=NULL))
      
      #download.file(paste0("https://drive.google.com/uc?id=",compoundMetaData$PDB_GID,"&export=download"), destfile = paste0("./data/",GENEID,'/PDB.zip'), mode = "wb")
      #Because PDB folder is very large, skipping virus warning
      ####gdown$download(paste0("https://drive.google.com/uc?id=",compoundMetaData$PDB_GID,"&export=download"), output = paste0("./data/",GENEID,'/PDB.zip'))
      unzip("KinCo_sel_2022-08-14.zip", files = c(paste0("KinCo_sel_2022-08-14/prot_PDB/", GENEID,".zip")), exdir = "data")
      pdb_needed <- unique(posedf$ModelName)
      unzip(paste0("data/KinCo_sel_2022-08-14/prot_PDB/", GENEID,".zip"), files = unlist(
        lapply(pdb_needed,function(m){paste0(m,".pdb.zip")})
      ), exdir = "data")
      #Deep pulling zipped data
      lapply(pdb_needed,function(m){
        zf <- zipfile$ZipFile(paste0("data/",m,".pdb.zip"),"r")
        zf$extract(unzip(paste0("data/",m,".pdb.zip"), list=T)$Name, path = "./data/")
        file.move(paste0("data/",unzip(paste0("data/",m,".pdb.zip"), list=T)$Name),"./data")
        unlink("data/n", recursive = TRUE)
      })
      
      showModal(modalDialog("Compound found!", footer=NULL))
      showModal(modalDialog(HTML("Locating our database... <br>
                            Kinase found! <br>
                            Searching for compound... <br>
                            Compound found! <br>
                            Loading the pair into server..."), footer=NULL))
      ###unzip(paste0("./data/",GENEID,'/PDB.zip') ,exdir=paste0("./data/",GENEID,'/'))
      #showModal(modalDialog(as.character(list.files(paste0("./data/",kinase_compound_pair[1],'/PDB/'))[1])))
      
      
      
    
      
      
      
      removeModal()
      
      
      
      
      observeEvent(input$addItem, {
        
        #if (!dir.exists(paste0("./data/",kinase_compound_pair[1],'/PDB/'))) {showNotification("Please wait ... Model still loading ...", type = "error")} 
        #else{
        counter$count=counter$count+1
        inserted <<- c(counter$count,inserted)
        if (counter$count %%3 ==1){
          insertUI(selector=paste0("#",ns("homologModulePlaceholder_first")),where="beforeEnd",
                   ui = tags$div(id = counter$count,
                                 homologInnerUi(id=ns(paste0("homologModule", counter$count ))))
          )}else if(counter$count %%3 ==2){
          insertUI(selector=paste0("#",ns("homologModulePlaceholder_second")),where="beforeEnd",
                   ui = tags$div(id = counter$count,
                                 homologInnerUi(id=ns(paste0("homologModule", counter$count ))))
          )}else{
          insertUI(selector=paste0("#",ns("homologModulePlaceholder_third")),where="beforeEnd",
                   ui = tags$div(id = counter$count,
                                 homologInnerUi(id=ns(paste0("homologModule", counter$count ))))
          )}
        
        homologServer(id=paste0("homologModule", counter$count ))
        
      }
      #}
      )
      
      observeEvent(input$removeItem, {
        
        removeUI(selector = paste0("#",inserted[length(inserted)]))
        inserted <<- inserted[-length(inserted)]
      }
      )
      
      
    }
    
  )
}
#####server code homolog UI

homologInnerUi<-function(id){
  
  #split data according to model names
  model_pairs <<- split.data.frame(posedf, f = posedf$ModelName)
  model_names <- names(model_pairs)
  
  
  ns=NS(id)
  
  fluidRow(
    
    bsCollapse(id = "collapsed_panel", open = "Homolog view", multiple = FALSE,
               bsCollapsePanel("Homolog view", style = "primary",
                               
    pickerInput(  inputId=ns("homologSelector"),
                  label = "Select Homolog",
                  choices= model_names,
                  selected = NULL,
                  multiple = FALSE 
                ),
    selectizeInput(  inputId=ns("poseSelector"),
                     label = "Select Multiple Pose(s)",
                     choices= NULL,
                     selected = NULL,
                     multiple = TRUE,
                     options = list(
                       plugins = list("remove_button")
                     )
    )
               ),

    bsCollapsePanel("Edit structure view",
    
    colourpicker::colourInput(
      inputId = ns("set_background_color"),
      label = "Set background color",
      closeOnClick = TRUE,
      value = "#FFFFFF"
    ),
    actionButton(
      inputId = ns("zoom_in"),
      label = "Zoom in",
      icon = icon("plus")
    ),
    actionButton(
      inputId = ns("zoom_out"),
      label = "Zoom out",
      icon = icon("minus")
    ),
    actionButton(
      inputId = ns("spin"),
      label = "Spin",
      icon = icon("sync-alt")
    ),
    actionButton(
      inputId = ns("unspin"),
      label = "Unspin",
      icon = icon("stop")
    ),
    actionButton(
      inputId = ns("clear"),
      label = "Clear",
      icon = icon("trash-alt")
    ),
    selectInput(
      inputId = ns("set_style"),
      label = "Set Style",
      choices = c("Stick","Line", "Cross", "Sphere", "Cartoon"),
      selected = "Cartoon"
    ),
    sliderInput(
      inputId = ns("set_slab"),
      label = "Set slab of view",
      min = -150,
      value = c(-50, 50),
      animate = TRUE,
      step = 5,
      max = 150,
      dragRange = TRUE
    ),
    radioButtons(
      inputId = ns("set_projection"),
      label = "Set view projection scheme",
      choices = c("perspective", "orthographic"),
      inline = TRUE
    ),
    sliderInput(
      inputId = ns("set_perceived_distance"),
      label = "Set perceived distance",
      min = 0,
      max = 500,
      value = 300
    )
    )),
    br(),
    r3dmolOutput(outputId = ns("model"), height = "500px", width = "100%")
  )
  
  
  
  
}

#updates
homologServer<-function(id){
  
  
  render_homolog<- function(homolog, poses){ 
    # Here input the homolog .pdb filepath, pose .pdb file path
    docked_model <- bio3d_docked(homolog, poses)
    #                             lapply(pose, function(x) paste("data", x, sep="/")))
    renderR3dmol({
      return(
        r3dmol(
          cartoonQuality = 20,
          lowerZoomLimit = 50,
          upperZoomLimit = 350,
          backgroundColor = "#FFFFFF"
        ) %>%
          m_add_model(data = m_bio3d(docked_model)
                      , format = "pdb") %>%
          #m_center() %>%
          m_zoom_to() %>%
          m_set_style( 
            style = m_style_cartoon( color = "#839DC7")) %>%
          m_set_style(
            sel = m_sel(ss = "s"),
            style = m_style_cartoon(color =  "#839DC7", arrows = TRUE)
          ) %>%
          m_set_style(
            sel = m_sel(resn = "UNL"),
            style = m_style_stick()
          ) %>%
          m_set_style(
            sel = m_sel(resn = "LIG"),
            style = m_style_stick()
          ) 
      )
    })
  }
  
  
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      
      #observe({
        
        # updatePickerInput(session,
        #                   inputId = "homologSelector",
        #                   choices = unique(posedf$ModelName)
        # )
      #})
      
      observeEvent(eventExpr = input$homologSelector, handlerExpr = {
      updateSelectizeInput(
          session,
          inputId="poseSelector",
          choices = model_pairs[[input$homologSelector]]$PoseID
        )
        
      })
      
      observeEvent(input$poseSelector, {
        
        session_ns <- session$ns('tmp')
        mod_id <- substr(session_ns, 1, nchar(session_ns)-4)
        
        ###pdbqt <- paste0(pairdir,'/PDBQT/',lgdinchikey,'.pdbqt') #read pose coordinate
        pdbqt <- paste0(pairpath,'/ligand.pdbqt')
        if (!dir.exists(paste0(output_dir,tgtgeneid,'-',lgdinchikey))){
          dir.create(paste0(output_dir, tgtgeneid, '-', lgdinchikey))
        }
        model_name<- input$homologSelector
        pose_ids <- input$poseSelector
        for (pose_id in pose_ids){ 
        
        ligandchain <- posedf %>% filter(ModelName == model_name, PoseID == pose_id) %>%.$LigandChain%>%.[[1]]
        ligandchain <- data.frame(x = ligandchain$x, y = ligandchain$y, z = ligandchain$z)
        posepath <- paste0(output_dir, tgtgeneid, '-', lgdinchikey, '/',model_name,'-',pose_id,'_lgd.pdb')
        if (!dir.exists(posepath)){
        export_PDBQT(pdbqt, ligandchain, posepath)}
        }
        #selected_poses<-input$poseSelector
        poses<- lapply(pose_ids, function(POSE) paste0(output_dir, tgtgeneid, '-', lgdinchikey, '/',model_name,'-',POSE,'_lgd.pdb'))
        ###homolog <- paste0("./data/",GENEID,'/PDB/',input$homologSelector,".pdb")
        homolog <- paste0("./data/",input$homologSelector,".pdb")
        
        output$model <- render_homolog(homolog, poses )  
      })
      
      observeEvent(input$set_background_color,{
        m_set_background_color(
          id = "model",
          hex = input$set_background_color
        )
      })
      
      observeEvent(input$spin, {
        m_spin(id = "model")
      })
      
      observeEvent(input$unspin, {
        m_spin("model", speed = 0)
      })
      
      observeEvent(input$zoom_out, {
        
        m_zoom(
          id = "model",
          factor = 0.7,
          animationDuration = 100)
        
      }) 
      
      observeEvent(input$zoom_in, {
        
        m_zoom(
          id = "model",
          factor = 1.3,
          animationDuration = 100)
        
      })
      
      observeEvent(input$set_projection, {
        
        m_set_projection(id = "model", scheme = input$set_projection)
      })
      
      observeEvent(input$clear, {
        m_clear(id = "model")
      })
      
      observeEvent(input$set_slab, {
        m_set_slab(
          id = "model",
          near = input$set_slab[1],
          far = input$set_slab[2]
        )
      })
      
      observeEvent(input$set_perceived_distance, {
        m_set_preceived_distance(id = "model", dist = input$set_perceived_distance)
      })
      
      observeEvent(input$set_style, {
        style <- switch(
          input$set_style,
          "Line" = list(line = list()),
          "Cartoon" = "Cartoon",
          "Stick" = list(stick = list()),
          "Cross" = list(cross = list()),
          "Sphere" = list(sphere = list())
        )
        
        if (input$set_style == "Cartoon"){
          
          m_set_style( id = "model",
            style = m_style_cartoon( color = "#839DC7"))
          
            m_set_style( id = "model",
              sel = m_sel(ss = "s"),
              style = m_style_cartoon(color =  "#839DC7", arrows = TRUE))
        
          
          
          m_set_style(
            id = "model",
            sel = m_sel(resn = "LIG"),
            style = m_style_stick()
          )
          m_set_style(
            id = "model",
            sel = m_sel(resn = "UNL"),
            style = m_style_stick()
          )
        } else{
        m_set_style(
          id = "model",
          style = style
        ) 
          m_set_style(
            id = "model",
            sel = m_sel(resn = "UNL"),
            style = m_style_stick()
          )
          m_set_style(
            id = "model",
            sel = m_sel(resn = "LIG"),
            style = m_style_stick()
          )
          } 
      })
    }
  )
}





##########server code - outer UI
# Initialize empty vector
inserted<- c()

#mainServer<-  function(id,data){
mainServer<-  function(id){
  moduleServer(
    id,
    function(input, output, session) {
      
      
      counter<-reactiveValues()
      
      counter$count=0
      
      ns <-session$ns
      
      
      observeEvent(input$kinaseSelector,{
        updateSelectizeInput(
          session,
          inputId="compoundSelector",
          choices = unique(kinco_meta[kinco_meta$Symbol==input$kinaseSelector,]$INCHIKEY_DESALT) %>% str_sub(1,-5)


        )
       
      })
      
      
      observeEvent(input$compoundSelector,{
        session_ns <- session$ns('tmp')
        mod_id <- substr(session_ns, 1, nchar(session_ns)-4)
      output$affinity <-renderText({kinco_meta[(kinco_meta$Symbol==input$kinaseSelector)&(kinco_meta$INCHIKEY_DESALT == paste0(input$compoundSelector,".zip")),]$affinity})
      
      
      })
      
      
      
      
      observeEvent(input$submit,{
        
        session_ns <- session$ns('tmp')
        mod_id <- substr(session_ns, 1, nchar(session_ns)-4)
        if (input$kinaseSelector =="") {showNotification("Please select a kinase to start with!", type = "error")} 
        else{removeUI(selector=paste0('#',ns('submit')), immediate=TRUE)
          insertUI(selector=paste0("#",ns("kinasePlaceholder")),where="beforeEnd",
                   ui = tags$div(
                     kinaseUi(id=ns(paste0("kinaseModule")))))
          
          
          #kinaseServer(id=paste0("kinaseModule", counter$count), data )
          kinaseServer(id=paste0("kinaseModule") , list(input$kinaseSelector, input$compoundSelector))
          }
      
        
        
      })
      
      
      #return(list(input$kinaseSelector, input$compoundSelector))
      
    }
    
  )
}




#mainUI

ui <- fluidPage(
  uiOutput("Module")
)

# main server
server <- function(input, output, session) {
  
  
  #data <- reactive({
    #return(data.frame(groups))
    #return(groups)
  #})
  
  output$Module <-renderUI({
    mainUI(id="firstTime" ) 
    
  })
  #mainServer(id="firstTime", data() )
  mainServer(id = "firstTime")
}

# run app
shinyApp(ui, server)
