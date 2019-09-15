# supportFunctions.R

source('~/code/R/gariUtil.R')

# Read the excel files with the tvalues and give back a processed data.table
myReadTvalues <- function(p) {


  glmNames           = p[['glmNames']]
  basedir            = p[['basedir']]
  vGuardarExcel      = p[['vGuardarExcel']]
  anArea             = p[['anArea']]
  bigSmall           = p[['bigSmall']]
  useFsaverageMNI152 = p[['useFsaverageMNI152']]
  Tmin               = p[['Tmin']]
  sujElimArtRepair   = p[['sujElimArtRepair']]

  letter = 'a'
  calculatePlots = F
  writeSvgPlots = F
  writeLatexPlots = F
  Obtain_BASIC_STATS = F
  writeExcel = F


  # The MNI coordinates can be calculated from fsaverage or from indiv space.
  # First I calculated from individual, but for comparison and to overlay labels
  # fsaverage is better. I have asked in the freesurfer list which one is more accurate
  # Meanwhile I am using fsaverage (305).
  if (useFsaverageMNI152){
    xGrep = '_fsx$'
    yGrep = '_fsy$'
    zGrep = '_fsz$'
  } else {
    xGrep = '_x$'
    yGrep = '_y$'
    zGrep = '_z$'
  }

  resultados = list()
  # glmName = 'block'
  for (glmName in glmNames){
    # cat('\n\n ***********************', glmName, '******************\n')

    DATA = list()

    sheetName = glmName
    baseMINIimgDir = file.path(basedir, 'imagesFromR')
    MINIimgDir = file.path(baseMINIimgDir,
                           paste0(anArea, glmName, vGuardarExcel,'Tmin',Tmin,letter))
    # dir.create(MINIimgDir, showWarnings = FALSE)

    # For reading:
    folder = paste(glmName,anArea, vGuardarExcel, sep="_")
    #??Listar todos los csv en el folder
    csvdir = file.path(csvBaseDir, folder)
    Datos_Entrada = dir(csvdir, pattern = '*.csv')
    DF = data.frame()
    for (di in Datos_Entrada){
      # The excel files where created with the function myGetLocalGlobalMaxima.m
      # This function reads the output files from the SPM analysis and creates excel files for further analysis
      data = read.csv(file.path(csvdir,di), header = T, sep = ',', row.names = 1,
                       dec = '.', fill = TRUE, numerals = 'allow.loss', na.strings = 'NaN')
      data = data[,-grep('X', names(data))]
      # for (ii in c(1:length(data))){
      #   data[,ii] = gsub('NaN', 'NA', data[,ii])
      #   data[,ii] = as.numeric(data[,ii])
      # }
      con <- gsub('(.*)_(.*)_(.*)_(.*)_(.*)', "\\2",di)
      con <- gsub('-', ".",con)
      data$Contrast = rep(con, nrow(data))
      data$SUBJECT = rownames(data)
      DF <- rbind(DF, data)
    }


    # It seems that subject S067 it is not already in the excel files, but to be sure do this:
    DF = DF[DF$SUBJECT %!in% "S067", ]
    # Then assigns metadata to each subjet
    GROUP = c(rep('ONCE', 35), rep('TWICE1', 29),
              rep('TWICE2', 30), 'TWICE1', 'TWICE2')
    PAIR  = c(rep('SOLO', 35),
              paste0('PAIR0', c(1:9)), paste0('PAIR', c(10:26,28:30)),
              paste0('PAIR0', c(1:9)), paste0('PAIR', c(10:30)),
              'PAIR31', 'PAIR31')
    TRT = c(rep('TEST', 64), rep('RETEST', 30), 'TEST', 'RETEST')
    DF$GROUP = as.factor(GROUP)
    DF$PAIR = as.factor(PAIR)
    DF$TRT = as.factor(TRT)

    #Remove the subjects identified as wrong
    DFS = DF[DF$SUBJECT %!in% sujElimArtRepair[[glmName]], ]
    # Select only the data that interests us
    selectInfo = names(DFS)[grep(paste(paste0(bigSmall,xGrep,'$'),
                                       paste0(bigSmall,yGrep,'$'),
                                       paste0(bigSmall,zGrep,'$'),
                                       paste0(bigSmall,'_T$'),
                                       '_vtx$', '_vtx305$',
                                       'Contrast',
                                       'SUBJECT','GROUP','PAIR','TRT',
                                       'inParc$','inVof$', sep='|'), names(DFS))]

    # Select only the contrasts we will report in the paper
    selectCon = c(#??'CSvsFF',    # Thesen 2012 (para letters), Price 1996
      # 'PWvsCS',    # James 2005
      # 'PWvsNull',  # Nestor 2013
      #??'PWvsRW',    # Haagoort 1999, Paulesu 2000, Xu 2001
      'RWvsCB',    # Bouhali 2014, Cohen 2002
      'RWvsCS',    # Thesen 2012, Cohen 2002-06, Rees 1999
      'RWvsFF',    # Woodhead 2011, Price 1996
      'RWvsPS',    # Ben-Shachar 2007-11, Rauschecker 2012,
      # Wang 14, Yeatman 2013-14
      'RWvsPW',    # Boukrina 2013, Graves 2010,
      # Cappa 1998, Kronbichler 2004
      'RWvsSD',    # Szwed 2011 > decir eliminamos
      # res parecidos PS, peor pq + esquinas
      # 'CSvsCB',    # Cohen 2002
      # 'RWvsFaces',  # Gauthier 2000
      # "FacesvsPSFaces",  # para ubicar el resto al lado de la faces.
      'RWvsNull')  # Binder 2006, Cohen 2000-08, Dehaene 2001-02-10,
    # Duncan 2009, Gelezer 2009-13-15, Longcamp 2011,
    # Nestor 2013, Fiez and Balota 1999,
    # Twomey 2011, Wang 2014, Wright 2008



    # Filter the dataset
    miniDF = DFS[DFS$Contrast %in% selectCon, ]
    miniDF = miniDF[,selectInfo]

    # Make a copy (I think there was code before this that I deleted, it is not useful)
    DFconNAth = copy(miniDF)


    # IF the T is not valid, then make the T value and the x,y,z values NaN so that they don't go to the means
    for (i in names(DFconNAth)[grep('_T', names(DFconNAth))]){
      outliers = !DFconNAth[, i] < Tmin
      outliers[outliers==FALSE] = NA
      DFconNAth[, i] = DFconNAth[, i] * outliers
      ROI = gsub('(.*)_(.*)', "\\1",i)
      if (useFsaverageMNI152){
        j = paste(ROI, 'fsx', sep='_')
        k = paste(ROI, 'fsy', sep='_')
        l = paste(ROI, 'fsz', sep='_')
      }else{
        j = paste(ROI, 'x', sep='_')
        k = paste(ROI, 'y', sep='_')
        l = paste(ROI, 'z', sep='_')
      }
      m = paste(ROI, 'vtx', sep='_')
      n = paste(ROI, 'vtx305', sep='_')
      DFconNAth[, k] = DFconNAth[, k] * outliers
      DFconNAth[, j] = DFconNAth[, j] * outliers
      DFconNAth[, l] = DFconNAth[, l] * outliers
      DFconNAth[, m] = DFconNAth[, m] * outliers
      DFconNAth[, n] = DFconNAth[, n] * outliers
    }



    # ORDER OF CONTRASTS AND TYPE
    DFconNAth$Contrast = factor(DFconNAth$Contrast)
    tiposContraste = c('Perceptual', "Lexical", "Lexical", 'Null',
                       "Perceptual", 'Lexical','Perceptual')

    DFconNAth$TYPE = NA
    for (ii in 1:length(tiposContraste)){
      Contraste = levels(DFconNAth$Contrast)[ii]
      Tipo      = tiposContraste[ii]
      DFconNAth[DFconNAth$Contrast == Contraste, c('TYPE')] = Tipo
    }
    # Make TYPE and subject a factor
    DFconNAth$TYPE = as.factor(DFconNAth$TYPE)
    DFconNAth$SUBJECT = as.factor(DFconNAth$SUBJECT)
    # Make the contrast the first column
    DFconNAth = DFconNAth[,c(grep('Contrast', names(DFconNAth)),
                             grep('Contrast', names(DFconNAth), invert = T))]


    if (0) { # FOLD THIS
      if (calculatePlots) {
        # plotHist    = list()
        plotScat    = list()
        # plotSubject = list()
        plotPair    = list()
        # plotBoxplot = list()
        # Create the plots
        #??Reorder the factors

        for (datoY in myGrepNames('^o', DFconNAth, T)[grep(
          yGrep, myGrepNames('^o', DFconNAth, T))]) {
          datoT = gsub(yGrep, '_T', datoY)
          datoX = gsub(yGrep, substr(xGrep,1,4), datoY)
          datoZ = gsub(yGrep, substr(zGrep,1,4), datoY)
          dato  = gsub(yGrep, '', datoY)

          # Para VOT usamos X e Y
          # Para IFG usamos Y y Z
          # Para PPC usamos X y Z
          USAR = c (datoZ, datoY)
          #??Scatterplot
          #chull function does not work with missing data
          nomissing <- na.omit(DFconNAth[,
                                         c('Contrast', 'TRT', 'GROUP', 'PAIR',
                                           'SUBJECT', USAR[1], USAR[2], datoT)])
          nomissing3 <- na.omit(DFconNAth[,
                                          c('Contrast', 'TRT', 'GROUP', 'PAIR',
                                            'SUBJECT', datoX, datoY, datoZ,datoT)])
          #getting the convex hull of each unique point set
          df <- nomissing
          find_hull <- function(df) df[chull(df[, USAR[1]], df[,USAR[2]]), ]
          hulls <- ddply(df, .(Contrast, TRT), find_hull)

          cat('plotScat', glmName, USAR[1], USAR[2],'...')

          plotScat[[USAR[1]]] <-
            ggplot(data=nomissing,
                   aes(x=get(USAR[1]),y=get(USAR[2]),colour=Contrast,
                       fill=Contrast, size=get(datoT))) +
            geom_point(data=labelPoints[labelPoints$ROI==dato,],
                       aes(x=x,y=y,colour=ROI,fill=ROI),alpha=0.2,size=.5) +
            geom_point() +
            geom_boxplot(aes(y=get(USAR[2]),  alpha=0.2)) +
            geom_boxplot(aes(y=get(USAR[1]),  alpha=0.2)) +
            geom_polygon(data = hulls, alpha = 0.25) +
            # geom_label(aes(label=SUBJECT)) +
            xlab(USAR[1]) + ylab(USAR[2]) + ggtitle(paste0('Tmin=',Tmin)) +
            facet_grid(TRT ~ Contrast) +
            theme_bw() +
            theme(text  =element_text(size=10),
                  axis.text=element_text(size=10),
                  # legend.position=c(.90, .90),
                  axis.title=element_text(size=12,face="bold"),
                  title=element_text(size=14, face="bold"),
                  strip.text = element_text(size = 12, angle = 0,
                                            face="bold"))
          cat('... done \n')


          plot_ly(lhIFG152, x = ~fsx, y = ~fsy, z = ~fsz,
                  marker = list(color = ~fsz,
                                colorscale = c('#FFE1A1', '#683531'),
                                showscale = TRUE)) %>% add_markers()
          subplots = list()
          for (Con in levels(nomissing3$Contrast)){
            subplots[[Con]] <- plot_ly(nomissing3[nomissing3$Contrast==Con,],
                                       x = ~lhIFG16_fsx,
                                       y = ~lhIFG16_fsy,
                                       z = ~lhIFG16_fsz,
                                       marker = list(color = ~lhIFG16_T,
                                                     colorscale = c('#FFE1A1', '#683531'),
                                                     showscale = TRUE)) %>%
              add_markers() %>%
              layout(scene = list(xaxis = list(title = 'X'),
                                  yaxis = list(title = 'Y'),
                                  zaxis = list(title = paste0(Con,'-Z'))))
          }
          subplot(subplots[["RWvsPS"]], subplots[["RWvsSD"]], subplots[["RWvsCB"]],
                  subplots[["RWvsCS"]], subplots[["RWvsFF"]], subplots[["RWvsPW"]],
                  nrows = 2)



          #??PAIR: Quiero ver solo TWICE1 y 2 y con linea entre mismo sujeto
          # por contraste
          solotwice = DFconNAth[DFconNAth$GROUP != 'ONCE', myGrepNames(
            '^out', DFconNAth, T)]
          logico = solotwice[is.na(solotwice[,USAR[2]]),c('Contrast', 'PAIR')]
          if (nrow(logico) != 0){
            for (ii in 1:nrow(logico)){
              solotwice[solotwice$Contrast==logico[ii,1] &
                          solotwice$PAIR==logico[ii,2], USAR[2]] = NA
            }
          }
          cat('plotPair', glmName, USAR[1], USAR[2],'...')
          plotPair[[USAR[2]]] =
            ggplot(data=solotwice,aes(x=PAIR,y=get(USAR[2]),colour=Contrast,
                                      label=SUBJECT))+
            geom_line(aes(group=PAIR)) + #, position=position_dodge(2)) +
            geom_point() +
            # geom_label(aes(size=get(datoT)/4)) +
            xlab(USAR[2]) +??ylab('SUBJECT') +
            ggtitle (paste0('TEST-RETEST anterior-posterior VOT,
                            Tmin=',Tmin)) +
            facet_grid(Contrast ~ ., scales = 'free_x') +
            theme_bw() +
            theme(text       =element_text(size=8),
                  axis.text  =element_text(size=5),
                  axis.text.x=element_text(angle=90),
                  axis.title =element_text(size=10,face="bold"),
                  title      =element_text(size=15, face="bold"),
                  strip.text =element_text(size=10,angle=0,face="bold"))
          cat('... done \n')
        }  # Fin FOR create plots
      }  # Fin IF calculatePlots??
      # Print the plots
      if (writeSvgPlots){
        width  = 20
        height = 20
        PS     = 12
        for (datoY in myGrepNames('^o', DFconNAth, T)[grep(
          yGrep, myGrepNames('^o', DFconNAth, T))]){

          datoT = gsub(yGrep, '_T', datoY)
          datoX = gsub(yGrep, xGrep, datoY)
          datoZ = gsub(yGrep, zGrep, datoY)
          dato  = gsub(yGrep, '', datoY)
          # ScatterplotUSAR
          cat('Write plotScat', glmName, USAR[1], USAR[2],'...')
          svg(filename = file.path(MINIimgDir,
                                   paste0("GROUPScat_", USAR[2],"_",
                                          vGuardarExcel, ".svg")),
              width = 18, height = 12, pointsize = 12)
          print(plotScat[[USAR[2]]])
          dev.off()
          cat('... done \n')

          # Test retest
          cat('Write plotPair', glmName, USAR[1], USAR[2],'...')
          svg(filename = file.path(MINIimgDir,
                                   paste0("GROUPPairs_", USAR[2],"_",
                                          vGuardarExcel, ".svg")),
              width = 18, height = 12, pointsize = 10)
          print(plotPair[[USAR[2]]])
          dev.off()
          cat('... done \n')
        }
      }
      if (writeLatexPlots) {
        for (datoY in myGrepNames('^o', DFconNAth, T)[grep(
          yGrep, myGrepNames('^o', DFconNAth, T))]){
          datoT = gsub(yGrep, '_T', datoY)
          USAR[1] = gsub(yGrep, 'x', datoY)
          # Scatterplot
          cat('Write plotScat', glmName, USAR[1], USAR[2],'...')
          print(plotScat[[USAR[1]]])
          cat('... done \n')

          # Test retest
          cat('Write plotPair', glmName, USAR[1], USAR[2],'...')
          print(plotPair[[USAR[2]]])
          cat('... done \n')
        }
      }

      # Obtain BASIC STATS
      if (Obtain_BASIC_STATS) {
        STATS = list()
        DTconNA = as.data.table(DFconNAth)
        for (i in names(DTconNA)[grep('SUBJECT|Contrast|out|GROUP|PAIR|TRT|vtx',
                                      names(DTconNA), invert = T)]){
          TMP = copy(DTconNA)
          TMP[,SUBJECT:=NULL]
          OUT = TMP[,list(mean=mean(get(i), na.rm=T),
                          sd=sd(get(i), na.rm=T),
                          N = sum(!is.na(get(i)))),
                    by=Contrast]
          STATS[[i]] = as.data.frame(OUT)
        }
         # end STATS

        STDT = as.data.frame(STATS)
        # Quitar las columnas .contrast y .N pq son repetidas
        names(STDT)[c(1)] = c('CONTRASTE')
        STDT = STDT[, -grep('*.Contrast$', names(STDT))]
        # Nos quedamos con los contrastes mas interesantes
        # OJO, orden diferente que en Contrasts.py, aqui lo ha ordenado por abcdario
        selectContrast  = c(1:nrow(STDT))
        selectInfo = c(1:ncol(STDT))

        miniSTDT = STDT[selectContrast, selectInfo]

        #??WRITE EXCEL
        #??Create folders
        csvFolders = file.path(basedir, 'excelFromR', paste0('Tmin', Tmin,
                                                             '_', vGuardarExcel))
        dir.create(csvFolders)
        csvSubFolders = file.path(csvFolders, 'SUBJECTS')
        dir.create(csvSubFolders)
       }

      if (writeExcel){
        for (suj in DFconNAth$SUBJECT){
          write.table(DFconNAth[DFconNAth$SUBJECT==suj,
                                myGrepNames('^out|GROUP|SUBJECT|PAIR|TRT',
                                            DFconNAth, T)],
                      file.path(csvSubFolders,paste0(glmName,'_',suj,letter,'.csv')),
                      col.names=T, row.names=F, quote = F, append = F,
                      sep = ',', na='NaN')
        }


        write.xlsx(miniSTDT,
                   file.path(csvFolders,
                             paste0('acpc_lhFusiform_Results_Tmin',Tmin,
                                    '_', vGuardarExcel, letter,'.xlsx')),
                   sheetName=paste0(sheetName, '_AvgPerCon'),
                   col.names=TRUE, row.names=TRUE, append=TRUE, showNA=TRUE)
      }
    }  # FOLD THIS




    # Store results
    DFconNAth = DFconNAth[,myGrepNames('^out', DFconNAth, T)]
    DFconNAth$DESIGN = as.factor(rep(glmName, nrow(DFconNAth)))
    resultados[[glmName]] = DFconNAth
  }  # Close the glmName

  # JOIN DESIGNS, as said, there is no 'event'
  allWithMissing = rbind(resultados[['block']], resultados[['event']])

  # Make ROIs long in order to be able to compare later on
  longROI = copy(allWithMissing)
  # Remove vertex
  longROI = longROI[, myGrepNames('vtx', longROI, T)]


  # Flip column names so that they end in _ROI
  varying = myGrepNames(paste(yGrep,xGrep,zGrep,'_T$', sep='|'), longROI, F)
  varying = gsub('(.*)_(.*)', "\\2.\\1",varying)
  names(longROI)[grep(paste(yGrep,xGrep,zGrep,'_T$',sep='|'), names(longROI), F)] = varying
  # Convert to long
  longROI = reshape(longROI, varying = varying,
                    direction = 'long', sep = '.')
  # Remove column 'id'
  longROI = longROI[,myGrepNames('id', longROI, T)]
  names(longROI)[grep('time', names(longROI))] = c('ROI')
  # Make ROI a factor
  longROI$ROI = as.factor(longROI$ROI)
  # Separate Null from the rest of the contrasts (we will use Null only in supplementary materials)
  longROINull = longROI[longROI$Contrast=='RWvsNull',]
  longROI = longROI[longROI$Contrast!='RWvsNull',]
  # Remove the categories from factors that no longer exist
  longROINull$Contrast = droplevels(longROINull$Contrast)
  longROI$Contrast     = droplevels(longROI$Contrast)
  longROINull$TYPE     = droplevels(longROINull$TYPE)
  longROI$TYPE         = droplevels(longROI$TYPE)



return(list('longROI'     = longROI,
            'longROINull' = longROINull)
       )

}


# Create the average and clustering algorithms, and plot it, with and without tracts
myClustering <-function(p) {
    # Libraries
    library(ggrepel)
    # trt = 'TEST'

    noDT = p[['noDT']]
    trt  = p[['trt']]

    listToReturn = list()


    # Leer tractos
    crearProbMapDWI = T
    if(crearProbMapDWI) {
      # This data was generated in the server, and it is safe in my gDrive
      roidir = '~/gDrive/BCBL/PROYECTOS/MINI/ANALYSIS/ROIfsaverage'
      # Separate it from TEST or RETEST from the beginning with trt variable
      Files   = list.files(path = roidir, pattern = paste0('^avgTract.*TEST_.*_THRcon.*[VOF|PARC]\\.csv$'))


        # READ THE FILES
        dFileList = list()
      for (dFile in Files){
        dFileList[[dFile]] = read.csv(file.path(roidir, dFile),header = F)[,c(1,2,3,4)]
        names(dFileList[[dFile]]) = c('x','y','z','prob')
      }

      # Crear DF para luego poder plotear
      bothTEST   = cbind(
        rbind(dFileList[[Files[grep('parc16_305_TEST',Files)]]],
              dFileList[[Files[grep('vof16_305_TEST',Files)]]]),
        c(rep('parc',nrow(dFileList[[Files[grep('parc16_305_TEST',Files)]]])),
          rep('vof' ,nrow(dFileList[[Files[grep('vof16_305_TEST',Files)]]])))
      )
      bothTEST$TRT = 'TEST'
      names(bothTEST) = c('x', 'y', 'z', 'prob','TRACT', 'TRT')
      bothRETEST   = cbind(
        rbind(dFileList[[Files[grep('parc16_305_RETEST',Files)]]],
              dFileList[[Files[grep('vof16_305_RETEST',Files)]]]),
        c(rep('parc',nrow(dFileList[[Files[grep('parc16_305_RETEST',Files)]]])),
          rep('vof' ,nrow(dFileList[[Files[grep('vof16_305_RETEST',Files)]]])))
      )
      bothRETEST$TRT = 'RETEST'
      names(bothRETEST) = c('x', 'y', 'z', 'prob','TRACT', 'TRT')


      dfTractsinDesign = rbind(bothTEST, bothRETEST)
      dfTract = cbind(rbind(dfTractsinDesign, dfTractsinDesign),
                      c(rep('block',nrow(dfTractsinDesign)),rep('event',nrow(dfTractsinDesign))))
      names(dfTract) = c('x','y','z', 'prob','TRACT','TRT', 'DESIGN')
      dfTract$TRT = as.factor(dfTract$TRT)
      dfTract$DESIGN = as.factor(dfTract$DESIGN)
      # AHORA THRESHOLDEAMOS LOS VALORES
      # Por ahora no lo thresholdeamos:
      # dfTractTh = dfTract[dfTract$prob > mean(dfTract$prob),]
      dfTractTh = dfTract[dfTract$prob > 0,]

      # Creamos mapa de colores para ggplot:
      #         parc    vof   (min es 1 para todos)
      # TEST    23      30
      # RETEST  50      57
      # Por ahora lo hardcodeo pero seria facil automatizarlo con esto:
      # max(dfTract$prob[dfTract$TRACT=='vof' & dfTract$TRT=='RETEST'])
      NContrast = length(levels(longROI$Contrast))
      # NClusters = 2  # Viene de abajo, en howManyGroups = 2, tiene que concordar
      # parcMax = 50
      # vofMax = 57
      parcMax = max(dfTract$prob[dfTract$TRACT=='parc' & dfTract$TRT==trt])
      vofMax = max(dfTract$prob[dfTract$TRACT=='vof' & dfTract$TRT==trt])
      # De todas maneras da igual pq siempre sera menor que NContrast

      colorSample = brewer.pal(NContrast, "Set2")
      # Estos dos cogeran los valores 1,2,3,4,5,6 y 1,2 respectivamente
      # Ahora tenemos que hacer que los valores de los plots concuerden
      allRcolors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
      colorPARC1 <- colorRampPalette(c("yellow"   ,"orange"  ))(floor(parcMax/2))
      colorPARC2 <- colorRampPalette(c("orange"   ,"red4"    ))(ceiling(parcMax/2))
      colorVOF1  <- colorRampPalette(c("lightsteelblue1","blue"))(floor(vofMax/2))
      colorVOF2  <- colorRampPalette(c("blue"     ,"darkblue"))(ceiling(vofMax/2))
      plotColors = c(colorSample,
                     colorPARC1 , colorPARC2,
                     colorVOF1  , colorVOF2)

      plotValues = c(
        # First six duplicated with number and contrast names
        # Changes cluster 1 and 2 colors to match the hue of the tract
        # '1'="#66C2A5",'2'="#FC8D62",'3'="#8DA0CB",'4'="#E78AC3",'5'="#A6D854",'6'="#FFD92F",
        '1'="#66C2A5",'2'="tomato3",'3'="#8DA0CB",'4'="#E78AC3",'5'="#A6D854",'6'="#FFD92F",
        # Now the same colors with contrast names.
        'RWvsCB'="#66C2A5",'RWvsCS'="#FC8D62",'RWvsFF'="#8DA0CB",
        'RWvsPS'="#E78AC3",'RWvsPW'="#A6D854",'RWvsSD'="#FFD92F",
        # Now the same colors for the TYPE
        'Lexical'="#66C2A5",'Perceptual'="#FC8D62",
        # Now the same colors for the CLUSTER
        'Cluster1'="#66C2A5",'Cluster2'="#FC8D62",'Cluster3'="#8DA0CB",
        'Cluster4'="#E78AC3",'Cluster5'="#A6D854",'Cluster6'="#FFD92F",
        # From 7 to 56 is for parc: 50
        '7'="#FFFF00",'8'="#FFFB00",'9'="#FFF700",'10'="#FFF300",
        '11'="#FFF000",'12'="#FFEC00",'13'="#FFE800",'14'="#FFE400",'15'="#FFE100",
        '16'="#FFDD00",'17'="#FFD900",'18'="#FFD500",'19'="#FFD200",'20'="#FFCE00",
        '21'="#FFCA00",'22'="#FFC600",'23'="#FFC300",'24'="#FFBF00",'25'="#FFBB00",
        '26'="#FFB700",'27'="#FFB400",'28'="#FFB000",'29'="#FFAC00",'30'="#FFA800",
        '31'="#FFA500",'32'="#FFA500",'33'="#FA9E00",'34'="#F59700",'35'="#F09000",
        '36'="#EB8900",'37'="#E68200",'38'="#E27B00",'39'="#DD7400",'40'="#D86E00",
        '41'="#D36700",'42'="#CE6000",'43'="#C95900",'44'="#C45200",'45'="#C04B00",
        '46'="#BB4400",'47'="#B63D00",'48'="#B13700",'49'="#AC3000",'50'="#A82900",
        '51'="#A32200",'52'="#9E1B00",'53'="#991400",'54'="#940D00",'55'="#8F0600",
        '56'="#8B0000",
        # And now for vof
        '57'="#CAE1FF",'58'="#C2D8FF",'59'="#BBD0FF",'60'="#B3C8FF",
        '61'="#ACBFFF",'62'="#A4B7FF",'63'="#9DAFFF",'64'="#95A6FF",'65'="#8E9EFF",
        '66'="#8696FF",'67'="#7F8DFF",'68'="#7785FF",'69'="#707DFF",'70'="#6874FF",
        '71'="#616CFF",'72'="#5963FF",'73'="#525BFF",'74'="#4A53FF",'75'="#434BFF",
        '76'="#3B42FF",'77'="#343AFF",'78'="#2C32FF",'79'="#2529FF",'80'="#1D21FF",
        '81'="#1619FF",'82'="#0E10FF",'83'="#0708FF",'84'="#0000FF",'85'="#0000FF",
        '86'="#0000FA",'87'="#0000F6",'88'="#0000F2",'89'="#0000EE",'90'="#0000EA",
        '91'="#0000E6",'92'="#0000E2",'93'="#0000DD",'94'="#0000D9",'95'="#0000D5",
        '96'="#0000D1",'97'="#0000CD",'98'="#0000C9",'99'="#0000C4",'100'="#0000C0",
        '101'="#0000BC",'102'="#0000B8",'103'="#0000B4",'104'="#0000B0",'105'="#0000AC",
        '106'="#0000A8",'107'="#0000A3",'108'="#00009F",'109'="#00009B",'110'="#000097",
        '111'="#000093",'112'="#00008F",'113'="#00008B")
      # Se les podria poner nombres como aqui:
      # scale_colour_manual(name="",
      #                       values = c("1"="yellow", "2"="orange", "3"="red",
      #                                  "value3"="grey", "value2"="black"))
      # Pero me ha parecido que era mas trabajo, mejor mapear por numero y ya esta
      dfTractTh$probColor = dfTractTh$prob
      dfTractTh$probColor[dfTractTh$TRACT=='parc'] = dfTractTh$probColor[dfTractTh$TRACT=='parc'] +
        NContrast
      dfTractTh$probColor[dfTractTh$TRACT=='vof'] =  dfTractTh$probColor[dfTractTh$TRACT=='vof'] +
        NContrast +
        parcMax
      # Y por ultimo lo convertimos en factor
      dfTractTh$probColor = as.factor(dfTractTh$probColor)
      # Hasta aqui era para pintar vof y parc

    }

    # # Cambiar cada vez si queremos hacerlo en todo el mask o en los dos rois definidos por tract
    # # Esto no afecta a IFG, solo a pPC y vOT
    # noDT = copy(noDT_mask)
    # # noDT = copy(noDT_DWI)
    # # noDT = copy(noDTp)  # esto es con los tres ROIs
    # str(noDT)


    # TRT y ROI en el MAIN experiment seran 1, TRT = TEST y ROI = lhVOT16

    # By independent contrast
    noDT[ ,Tmean:=mean(T, na.rm=TRUE),    by=c('Contrast','DESIGN','TRT','ROI')]
    noDT[ ,xmean:=mean(fsx, na.rm=TRUE),  by=c('Contrast','DESIGN','TRT','ROI')]
    noDT[ ,ymean:=mean(fsy, na.rm=TRUE),  by=c('Contrast','DESIGN','TRT','ROI')]
    noDT[ ,zmean:=mean(fsz, na.rm=TRUE),  by=c('Contrast','DESIGN','TRT','ROI')]
    # The contrasts are aggregated by type
    noDT[ ,TTYPE:=mean(T, na.rm=TRUE),    by=c('TYPE',    'DESIGN','TRT','ROI')]
    noDT[ ,xTYPE:=mean(fsx, na.rm=TRUE),  by=c('TYPE',    'DESIGN','TRT','ROI')]
    noDT[ ,yTYPE:=mean(fsy, na.rm=TRUE),  by=c('TYPE',    'DESIGN','TRT','ROI')]
    noDT[ ,zTYPE:=mean(fsz, na.rm=TRUE),  by=c('TYPE',    'DESIGN','TRT','ROI')]
    # Some sd-s
    noDT[ ,Tsd  :=sd(T, na.rm=TRUE)  ,    by=c('Contrast','DESIGN','TRT','ROI')]
    noDT[ ,xsd  :=sd(fsx, na.rm=TRUE)  ,  by=c('Contrast','DESIGN','TRT','ROI')]
    noDT[ ,ysd  :=sd(fsy, na.rm=TRUE)  ,  by=c('Contrast','DESIGN','TRT','ROI')]
    noDT[ ,zsd  :=sd(fsz, na.rm=TRUE)  ,  by=c('Contrast','DESIGN','TRT','ROI')]
    noDT[ ,N    := (.N - sum(is.na(fsx))),by=c('Contrast','DESIGN','TRT','ROI')]
    noDT[ ,NTYPE:= (.N - sum(is.na(fsx))),by=c('TYPE'    ,'DESIGN','TRT','ROI')]


    # Save it to return it
    listToReturn[['noDT']] = noDT



    # Validar y recuperar la info
    # dft = unique(noDT[,list(DESIGN,TRT,Contrast,N,ROI,TYPE,xmean,ymean,zmean,Tmean,xsd,ysd,zsd,Tsd)])
    dft = unique(noDT[,list(DESIGN,TRT,Contrast,N,NTYPE,ROI,TYPE,
                            xmean,ymean,zmean,Tmean,
                            xTYPE,yTYPE,zTYPE,TTYPE,
                            xsd,ysd,zsd,Tsd)])
    dft$size      = (dft$xsd + dft$ysd)
    dft$LABEL     = paste0(dft$Contrast,'(N=',dft$N,')')
    dft$LABELTYPE = paste0(dft$TYPE,'(N=',dft$NTYPE,')')


    # CLUSTER ANALYSIS
    howManyGroups = 2
    dft[,CLUSTER := myAgrupaDT(xmean,ymean,zmean,Tmean,xsd,ysd,zsd,Tsd, howManyGroups,'y'),
        by=list(DESIGN, TRT,ROI)]
    dft$CLUSTER = as.factor(dft$CLUSTER)
    dft[,CLUSTERT := myAgrupaDT(xmean,ymean,zmean,Tmean,xsd,ysd,zsd,Tsd, howManyGroups,'T'),
        by=list(DESIGN, TRT,ROI)]
    dft$CLUSTERT = as.factor(dft$CLUSTERT)
    dft[,CLUSTERx := myAgrupaDT(xmean,ymean,zmean,Tmean,xsd,ysd,zsd,Tsd, howManyGroups,'x'),
        by=list(DESIGN, TRT,ROI)]
    dft$CLUSTERT = as.factor(dft$CLUSTERT)
    dft[,CLUSTERz := myAgrupaDT(xmean,ymean,zmean,Tmean,xsd,ysd,zsd,Tsd, howManyGroups,'z'),
        by=list(DESIGN, TRT,ROI)]
    dft$CLUSTERT = as.factor(dft$CLUSTERT)
    # print(unique(dft[,list(DESIGN, Contrast, xmean, ymean,zmean,Tmean,CLUSTER,CLUSTERT)]))


    # Save it to return it
    listToReturn[['dft']] = dft






    # PLOT
    lado = 13
    ver = 'VOTv06_paper'
    vOTplots = T
    if(vOTplots){
      versionSinTractos = T
      if(versionSinTractos){
        createPlotsXY = T
        if(createPlotsXY){
          # contrastPlotXY = ggplot( data=dft, aes(x=xmean,y=ymean)) +
          contrastPlotXY = ggplot( data=dft[dft$DESIGN=='block']) +
            # annotate('rect', xmin = aVWFA152[1]-(lado/2), xmax = aVWFA152[1]+(lado/2),
            #                  ymin = aVWFA152[2]-(lado/2), ymax = aVWFA152[2]+(lado/2),
            #                  fill = "blue", colour = 'grey', alpha=0.2) +
            # annotate('rect', xmin = cVWFA152[1]-(lado/2), xmax = cVWFA152[1]+(lado/2),
            #                  ymin = cVWFA152[2]-(lado/2), ymax = cVWFA152[2]+(lado/2),
            #                  fill = "blue", colour = 'grey', alpha=0.2) +
            # annotate('rect', xmin = pVWFA152[1]-(lado/2), xmax = pVWFA152[1]+(lado/2),
            #                  ymin = pVWFA152[2]-(lado/2), ymax = pVWFA152[2]+(lado/2),
            #                  fill = "blue", colour = 'grey', alpha=0.2) +
            annotate('text', y = aVWFA152[1]-1, x = aVWFA152[2],  # -(lado/3),
                     label='aVWFA',colour = 'gray25') +
            annotate('text', y = cVWFA152[1]-1, x = cVWFA152[2],  #??-(lado/3),
                     label='cVWFA',colour = 'gray25') +
            annotate('text', y = pVWFA152[1]-1, x = pVWFA152[2],  # -(lado/3),
                     label='pVWFA',colour = 'gray25') +
            annotate('text', y = aVWFA152[1], x = aVWFA152[2],  # -(lado/3),
                     label='+',colour = 'gray25') +
            annotate('text', y = cVWFA152[1], x = cVWFA152[2],  #??-(lado/3),
                     label='+',colour = 'gray25') +
            annotate('text', y = pVWFA152[1], x = pVWFA152[2],  # -(lado/3),
                     label='+',colour = 'gray25') +
            # geom_point(x = aVWFA152[1], y = aVWFA152[2],inherit.aes=F,
            #            shape=3, colour = 'gray25') +
            # geom_point(x = -55, y = -42, inherit.aes=F,shape=3, colour = 'gray25') +
            # geom_point(y = cVWFA152[1], x = cVWFA152[2],inherit.aes=F,
            #                shape=3, colour = 'gray25') +
            # geom_point(y = pVWFA152[1], x = pVWFA152[2],inherit.aes=F,
            #                shape=3, colour = 'gray25') +
            # geom_point(data = dfTractTh[dfTractTh$TRACT=='vof',],
            #             aes(x,y,color=probColor,fill=probColor),size=1.2) +
            # geom_point(data = dfTractTh[dfTractTh$TRACT=='parc',],
            #           aes(x,y,color=probColor,fill=probColor),size=1.2) +
          # scale_colour_manual(name="", values = plotValues) +
          scale_color_manual(name="", values = c("#0072B2",
                                                 "tomato3",
                                                 "Perceptual"="#0072B2",
                                                 "Lexical"="tomato3")) +
            geom_point(aes(y=xmean,x=ymean,size=size, fill=TYPE,colour=TYPE), alpha=0.7,inherit.aes=F) +

            scale_size_identity() +
            geom_point(aes(y=xmean,x=ymean,size=Tmean),
                       colour='black',fill='black', alpha=1, inherit.aes=F) +
            geom_label_repel(aes(y=xmean,x=ymean, label=LABEL),
                             size = 3, inherit.aes=F, colour='black') +

            xlab('MNI Y') + ylab('MNI X') +
            scale_x_reverse(breaks=seq(-75,-50,5), limits= c(-48,-75)) +
            # scale_x_continuous(breaks=seq(-12,-6,2), limits= c(-15,-3)) +
            # scale_y_continuous(breaks=seq(-80,-45,5), limits= c(-80,-44)) +
            scale_y_continuous(breaks=seq(-46,-40,2), limits= c(-47,-38)) +
            # ggtitle(paste0('Clustered Contrasts, Tmin: ', Tmin)) +
            # facet_grid(TRT ~ DESIGN) +
            # facet_grid(. ~ DESIGN) +
            facet_grid(TRT ~ .) +
            coord_equal() +
            theme_bw() +
            theme(text  =element_text(size=8),
                  legend.position="none",
                  axis.text=element_text(size=8),
                  # legend.position=c(.90, .90),
                  axis.title=element_text(size=10,face="bold"),
                  title=element_text(size=10, face="bold"),
                  strip.text = element_text(size = 10, angle = 0,face="bold"),
                  # panel.border = element_blank(),
                  # panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()
                  # axis.line = element_line(colour = "black")
            )
          # DO NOT PLOT IT, SAVE IT AND RETURN IT
          listToReturn[['contrastPlotXY']] = contrastPlotXY
          # Plot: print on screen, high q vectorial, low q png
          # print(contrastPlotXY)
          # svg(filename = file.path(PathFigsFldrCap7,'ContrastClusterAverages',
          #                          paste0("tractContrastPlotXY_",ver,".svg")),
          #     width = 5, height = 3, pointsize = 12)
          # # png(filename = file.path(PathFiguresFolder, paste0("contrastPlotXY_lq_",ver,".png")),
          # #                 width=400,height=550,units='px',pointsize=8,bg='white',res=90)
          # print(contrastPlotXY)
          # dev.off()
        }
        createPlotsZY = F
        if(createPlotsZY){
          contrastPlotZY = ggplot( data=dft, aes(x=zmean,y=ymean)) +
            # annotate('rect', xmin = aVWFA152[3]-(lado/2), xmax = aVWFA152[3]+(lado/2),
            #                  ymin = aVWFA152[2]-(lado/2), ymax = aVWFA152[2]+(lado/2),
            #                  fill = "blue", colour = 'grey', alpha=0.2) +
            # annotate('rect', xmin = cVWFA152[3]-(lado/2), xmax = cVWFA152[3]+(lado/2),
            #                  ymin = cVWFA152[2]-(lado/2), ymax = cVWFA152[2]+(lado/2),
            #                  fill = "blue", colour = 'grey', alpha=0.2) +
            # annotate('rect', xmin = pVWFA152[3]-(lado/2), xmax = pVWFA152[3]+(lado/2),
            #                  ymin = pVWFA152[2]-(lado/2), ymax = pVWFA152[2]+(lado/2),
            #                  fill = "blue", colour = 'grey', alpha=0.2) +
            annotate('text', x = aVWFA152[3]+(lado/3), y = aVWFA152[2]+(lado*0),
                     label='aVWFA',colour = 'gray25') +
            annotate('text', x = cVWFA152[3]+(lado/3), y = cVWFA152[2]+(lado*0),
                     label='cVWFA',colour = 'gray25') +
            annotate('text', x = pVWFA152[3]-(lado/3), y = pVWFA152[2]-(lado*0),
                     label='pVWFA',colour = 'gray25') +
            geom_point(x = aVWFA152[3], y = aVWFA152[2],inherit.aes=F,
                       shape=3, colour = 'gray25') +
            geom_point(x = cVWFA152[3], y = cVWFA152[2],inherit.aes=F,
                       shape=3, colour = 'gray25') +
            geom_point(x = pVWFA152[3], y = pVWFA152[2],inherit.aes=F,
                       shape=3, colour = 'gray25') +
            # geom_point(data = dfTractTh[dfTractTh$TRACT=='vof',],
            #             aes(x,y,color=probColor,fill=probColor),size=1.2) +
            # geom_point(data = dfTractTh[dfTractTh$TRACT=='parc',],
            #           aes(x,y,color=probColor,fill=probColor),size=1.2) +
            scale_colour_manual(name="", values = plotValues) +
            geom_point(aes(size=size, fill=CLUSTER,colour=CLUSTER), alpha=0.4) +
            scale_size_identity() +
            geom_point(aes(x=zmean,y=ymean,size=Tmean),
                       colour='black',fill='black', alpha=0.6,inherit.aes=F) +
            geom_label_repel(aes(x=zmean,y=ymean, label=LABEL),
                             size = 3, inherit.aes=F, colour='black') +
            xlab('MNI Z') + ylab('MNI Y') +
            scale_x_continuous(breaks=seq(-16,-4,2), limits= c(-17,-2)) +
            scale_y_continuous(breaks=seq(-80,-45,5), limits= c(-80,-44)) +
            # ggtitle(paste0('Clustered Contrasts, Tmin: ', Tmin)) +
            # facet_grid(TRT ~ DESIGN) +
            facet_grid(. ~ DESIGN) +
            coord_equal() +
            theme_bw() +
            theme(text  =element_text(size=8),
                  legend.position="none",
                  axis.text=element_text(size=8),
                  # legend.position=c(.90, .90),
                  axis.title=element_text(size=10,face="bold"),
                  title=element_text(size=10, face="bold"),
                  strip.text = element_text(size = 10, angle = 0,face="bold"),
                  # panel.border = element_blank(),
                  # panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()
                  # axis.line = element_line(colour = "black")
            )
          # DO NOT PLOT IT, SAVE IT AND RETURN IT
          # listToReturn[['contrastPlotZY']] = contrastPlotZY
          # # Plot: print on screen, high q vectorial, low q png
          # print(contrastPlotZY)
          # # print(contrastPlotZY)
          # png(filename = file.path(PathFiguresFolder, paste0("contrastPlotZY_lq_",ver,".png")),
          #     width=400,height=550,units='px',pointsize=8,bg='white',res=90)
          # print(contrastPlotZY)
          # dev.off()
        }
      } # Cierra versionSinTractos
      versionConTractos = T
      if(versionConTractos){
        createPlotsXY = T
        if(createPlotsXY){
          # contrastPlotXY = ggplot( data=dft, aes(x=xmean,y=ymean)) +

          contrastPlotXYcon = ggplot( data=dft[dft$DESIGN=='block',], aes(x=ymean,y=xmean)) +
            geom_point(data = dfTractTh[dfTractTh$TRACT=='vof',],
                       aes(y,x,color=probColor,fill=probColor),size=1.2) +
            geom_point(data = dfTractTh[dfTractTh$TRACT=='parc',],
                       aes(y,x,color=probColor,fill=probColor),size=1.2) +
            # annotate('text', y = aVWFA152[1]-1, x = aVWFA152[2],  # -(lado/3),
            #          label='aVWFA',colour = 'gray25') +
            # annotate('text', y = cVWFA152[1]-1, x = cVWFA152[2],  #??-(lado/3),
            #          label='cVWFA',colour = 'gray25') +
            # annotate('text', y = pVWFA152[1]-1, x = pVWFA152[2],  # -(lado/3),
            #          label='pVWFA',colour = 'gray25') +
            # annotate('text', y = aVWFA152[1], x = aVWFA152[2],  # -(lado/3),
            #          label='+',colour = 'gray25') +
            # annotate('text', y = cVWFA152[1], x = cVWFA152[2],  #??-(lado/3),
            #          label='+',colour = 'gray25') +
            # annotate('text', y = pVWFA152[1], x = pVWFA152[2],  # -(lado/3),
          #          label='+',colour = 'gray25') +
          scale_colour_manual(name="", values = plotValues) +
            geom_point(aes(size=size), shape=1, colour='gray25', alpha=0.4) +
            scale_size_identity() +
            geom_point(aes(x=ymean,y=xmean,size=Tmean),
                       colour='black',fill='black', alpha=0.6, inherit.aes=F) +
            geom_label_repel(aes(x=ymean,y=xmean, label=LABEL),
                             size = 3, inherit.aes=F, colour='black') +
            xlab('MNI Y') + ylab('MNI X') +
            scale_y_reverse(breaks=seq(-60,0,10), limits= c(0,-60)) +
            scale_x_reverse(breaks=seq(-100,-10,10), limits= c(-10,-102)) +
            # ggtitle(paste0('Clustered Contrasts, Tmin: ', Tmin)) +
            # facet_grid(TRT ~ DESIGN) +
            # facet_grid(. ~ DESIGN) +
            coord_equal() +
            theme_bw() +
            theme(text  =element_text(size=12),
                  legend.position="none",
                  axis.text=element_text(size=12),
                  # legend.position=c(.90, .90),
                  axis.title=element_text(size=14,face="bold"),
                  title=element_text(size=12, face="bold"),
                  strip.text = element_text(size = 16, angle = 0,face="bold"),
                  # panel.border = element_blank(),
                  # panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()
                  # axis.line = element_line(colour = "black")
            )
          # DO NOT PLOT IT, SAVE IT AND RETURN IT
          listToReturn[['contrastPlotXYcon']] = contrastPlotXYcon
          # Plot: print on screen, high q vectorial, low q png
          # print(contrastPlotXY)
          # svg(filename = file.path(PathFigsFldrCap7,'tractContrastClusterAverages',
          #                          paste0("tractContrastPlotXY_",ver,".svg")),
          #     width = 5, height = 3, pointsize = 12)
          # print(contrastPlotXY)
          # png(filename = file.path(PathFigsFldrCap7,'tractContrastClusterAverages',
          #                          paste0("tractContrastPlotXY_",ver,".png")),
          #     width=5,height=3,units='in',pointsize=12,bg='white',res=300)
          # print(contrastPlotXY)
          # dev.off()
        }
        createPlotsZY = F
        if(createPlotsZY){
          contrastPlotZY = ggplot( data=dft, aes(x=zmean,y=ymean)) +
            geom_point(data = dfTractTh[dfTractTh$TRACT=='vof',],
                       aes(z,y,color=probColor,fill=probColor),size=1.2) +
            geom_point(data = dfTractTh[dfTractTh$TRACT=='parc',],
                       aes(z,y,color=probColor,fill=probColor),size=1.2) +
            scale_colour_manual(name="", values = plotValues) +
            geom_point(aes(size=size, fill=CLUSTER,colour=CLUSTER), alpha=0.4) +
            scale_size_identity() +
            geom_point(aes(x=zmean,y=ymean,size=Tmean),
                       colour='black',fill='black', alpha=0.6,inherit.aes=F) +
            annotate('text', x = aVWFA152[3]-(lado/3), y = aVWFA152[2]+(lado/3),
                     label='aVWFA',colour = 'gray25') +
            annotate('text', x = cVWFA152[3]-(lado/3), y = cVWFA152[2]-(lado/3),
                     label='cVWFA',colour = 'gray25') +
            annotate('text', x = pVWFA152[3]-(lado/3), y = pVWFA152[2]-(lado/3),
                     label='pVWFA',colour = 'gray25') +
            geom_point(x = aVWFA152[3], y = aVWFA152[2],inherit.aes=F,
                       shape=3, colour = 'gray25') +
            geom_point(x = cVWFA152[3], y = cVWFA152[2],inherit.aes=F,
                       shape=3, colour = 'gray25') +
            geom_point(x = pVWFA152[3], y = pVWFA152[2],inherit.aes=F,
                       shape=3, colour = 'gray25') +
            geom_label_repel(aes(x=zmean,y=ymean, label=LABEL),
                             size = 3, inherit.aes=F, colour='black') +
            xlab('MNI Z') + ylab('MNI Y') +
            scale_x_continuous(breaks=seq(-30,30,10), limits= c(-30,30)) +
            scale_y_continuous(breaks=seq(-100,-10,10), limits= c(-102,-10)) +
            # ggtitle(paste0('Clustered Contrasts, Tmin: ', Tmin)) +
            # facet_grid(TRT ~ DESIGN) +
            facet_grid(. ~ DESIGN) +
            coord_equal() +
            theme_bw() +
            theme(text  =element_text(size=12),
                  legend.position="none",
                  axis.text=element_text(size=12),
                  # legend.position=c(.90, .90),
                  axis.title=element_text(size=14,face="bold"),
                  title=element_text(size=12, face="bold"),
                  strip.text = element_text(size = 16, angle = 0,face="bold"),
                  # panel.border = element_blank(),
                  # panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()
                  # axis.line = element_line(colour = "black")
            )
          # Plot: print on screen, high q vectorial, low q png
          # print(contrastPlotZY)
          # svg(filename = file.path(PathFigsFldrCap7,'tractContrastClusterAverages',
          #                          paste0("tractContrastPlotZY_",ver,".svg")),
          #     width = 6, height = 5, pointsize = 12)
          # print(contrastPlotZY)
          # png(filename = file.path(PathFiguresFolder, paste0("tractContrastPlotZY_lq_",ver,".png")),
          #                 width=600,height=750,units='px',pointsize=8,bg='white',res=90)
          # print(contrastPlotZY)
          # dev.off()
        }
      } # Cierra versionConTractos
    }

    # PLOTEAMOS CON LOS CONTRASTES AGRUPADOS EN TYPE/RESULTADO DEL CLUSTERING
    # Primero creamos el clustering basado en lo anterior
    clusterNames = c('Cluster1','Cluster2','Cluster3','Cluster4','Cluster5','Cluster6')
    noDFc = copy(noDT)
    noDFc$CLUSTER = 0
    # noDFc = as.data.frame(noDFc[noDFc$ROI==roi,])
    for (roi in levels(dft$ROI)){
      for (trt in levels(dft$TRT)){  # FOR TEST
        for (design in levels(dft$DESIGN)){
          matcheo = dft[dft$TRT==trt & dft$DESIGN==design & dft$ROI==roi,
                        list(Contrast,CLUSTER)]
          for (contrast in matcheo[,Contrast]){
            noDFc[noDFc$TRT==trt & noDFc$DESIGN==design & noDFc$Contrast==contrast,
                  'CLUSTER'] = as.double(matcheo[matcheo$Contrast==contrast, 'CLUSTER',with=F])
          }}}}
    noDFc$CLUSTER = as.factor(noDFc$CLUSTER)
    levels(noDFc$CLUSTER) = clusterNames[1:howManyGroups]

    noDTc = as.data.table(noDFc)
    noDTc[ ,cxmean:=mean(fsx, na.rm=TRUE), by=c('TYPE','DESIGN','TRT','ROI')]
    noDTc[ ,cymean:=mean(fsy, na.rm=TRUE), by=c('TYPE','DESIGN','TRT','ROI')]
    noDTc[ ,czmean:=mean(fsz, na.rm=TRUE), by=c('TYPE','DESIGN','TRT','ROI')]
    noDTc[ ,cTmean:=mean(T, na.rm=TRUE), by=c('TYPE','DESIGN','TRT','ROI')]
    noDTc[ ,cxsd  :=sd(fsx, na.rm=TRUE)  , by=c('TYPE','DESIGN','TRT','ROI')]
    noDTc[ ,cysd  :=sd(fsy, na.rm=TRUE)  , by=c('TYPE','DESIGN','TRT','ROI')]
    noDTc[ ,czsd  :=sd(fsz, na.rm=TRUE)  , by=c('TYPE','DESIGN','TRT','ROI')]
    noDTc[ ,cTsd  :=sd(T, na.rm=TRUE)  , by=c('TYPE','DESIGN','TRT','ROI')]
    noDTc[ ,cN    := .N                , by=c('TYPE','DESIGN','TRT','ROI')]

    df = unique(noDTc[,list(DESIGN,TRT,Contrast,cN,ROI,TYPE,cxmean,cymean,czmean,cTmean,
                            cxsd,cysd,czsd,cTsd,CLUSTER)])
    df$csize = (df$cxsd + df$cysd)
    df$LABEL = paste0(df$Contrast,'(',df$cN,')')

    # Save it to return it
    listToReturn[['dfc']] = df

    # print(df)
    # cat(roi, '\n')
    versionSinTractos = T
    if(versionSinTractos){
      createClusterPlotXY = T
      if(createClusterPlotXY){
        clusterPlotXY = ggplot(data=df[df$DESIGN=='block']) +
          annotate('text', y = aVWFA152[1]-1, x = aVWFA152[2],  # -(lado/3),
                   label='aVWFA',colour = 'gray25') +
          annotate('text', y = cVWFA152[1]-1, x = cVWFA152[2],  #??-(lado/3),
                   label='cVWFA',colour = 'gray25') +
          annotate('text', y = pVWFA152[1]-1, x = pVWFA152[2],  # -(lado/3),
                   label='pVWFA',colour = 'gray25') +
          annotate('text', y = aVWFA152[1], x = aVWFA152[2],  # -(lado/3),
                   label='+',colour = 'gray25') +
          annotate('text', y = cVWFA152[1], x = cVWFA152[2],  #??-(lado/3),
                   label='+',colour = 'gray25') +
          annotate('text', y = pVWFA152[1], x = pVWFA152[2],  # -(lado/3),
                   label='+',colour = 'gray25') +
          scale_color_manual(name="", values = c("#0072B2",
                                                 "tomato3",
                                                 "Perceptual"="#0072B2",
                                                 "Lexical"="tomato3")) +
          # scale_colour_manual(name="", values = plotValues) +
          geom_point(aes(x=cymean,y=cxmean, size=csize,colour=TYPE, fill=TYPE),
                     alpha=0.7,inherit.aes=F) +
          scale_size_identity() +
          # geom_label(aes(y=cxmean,x=cymean, label=TYPE),nudge_y = -2, nudge_x = 1,
          geom_label_repel(data=dft[dft$DESIGN=='block',], aes(x=yTYPE,y=xTYPE, label=LABELTYPE),
                     colour='black', size = 3, inherit.aes=F) +
          # geom_point(aes(size=cTmean),colour='black', fill='black', alpha = 0.7) +
          geom_point(aes(y=cxmean,x=cymean,size=cTmean),colour='black', fill='black', alpha = 0.7) +
          xlab('MNI Y') + ylab('MNI X') +
          scale_x_reverse(breaks=seq(-75,-50,5), limits= c(-48,-75)) +
          scale_y_continuous(breaks=seq(-46,-40,2), limits= c(-47,-38)) +
          # facet_grid(TRT ~ DESIGN) +
          facet_grid(TRT ~ .) +
          # facet_grid(. ~ DESIGN) +
          coord_equal() +
          theme_bw() +
          theme(text  =element_text(size=8),
                legend.position="none",
                axis.text=element_text(size=8),
                # legend.position=c(.90, .90),
                axis.title=element_text(size=10,face="bold"),
                title=element_text(size=10, face="bold"),
                strip.text = element_text(size = 10, angle = 0,face="bold"),
                # panel.border = element_blank(),
                # panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
                # axis.line = element_line(colour = "black")
          )
        # DO NOT PLOT IT, SAVE IT AND RETURN IT
        listToReturn[['clusterPlotXY']] = clusterPlotXY
        # # Plot: print on screen, high q vectorial, low q png
        # print(clusterPlotXY)
        # svg(filename = file.path(PathFigsFldrCap7,'ContrastClusterAverages',
        #                          paste0("CLUSTERContrastPlotXY_",ver,".svg")),
        #     width = 5, height = 3, pointsize = 12)
        # # png(filename = file.path(PathFiguresFolder, paste0("clusterPlotXY_lq_",ver,".png")),
        # #                 width=400,height=550,units='px',pointsize = 8, bg ='white',res=90)
        # print(clusterPlotXY)
        # dev.off()
      }
      createClusterPlotZY = F
      if(createClusterPlotZY){
        clusterPlotZY = ggplot(data=df,aes(x=czmean,y=cymean)) +
          # annotate('rect', xmin = aVWFA152[3]-(lado/2), xmax = aVWFA152[3]+(lado/2),
          #                  ymin = aVWFA152[2]-(lado/2), ymax = aVWFA152[2]+(lado/2),
          #                  fill = "blue", colour = 'grey', alpha=0.2) +
          # annotate('rect', xmin = cVWFA152[3]-(lado/2), xmax = cVWFA152[3]+(lado/2),
          #                  ymin = cVWFA152[2]-(lado/2), ymax = cVWFA152[2]+(lado/2),
          #                  fill = "blue", colour = 'grey', alpha=0.2) +
          # annotate('rect', xmin = pVWFA152[3]-(lado/2), xmax = pVWFA152[3]+(lado/2),
          #                  ymin = pVWFA152[2]-(lado/2), ymax = pVWFA152[2]+(lado/2),
          #                  fill = "blue", colour = 'grey', alpha=0.2) +
          annotate('text', x = aVWFA152[3]+(lado/3), y = aVWFA152[2]+(lado*0),
                   label='aVWFA',colour = 'gray25') +
          annotate('text', x = cVWFA152[3]+(lado/3), y = cVWFA152[2]+(lado*0),
                   label='cVWFA',colour = 'gray25') +
          annotate('text', x = pVWFA152[3]-(lado/3), y = pVWFA152[2]-(lado*0),
                   label='pVWFA',colour = 'gray25') +
          geom_point(x = aVWFA152[3], y = aVWFA152[2],inherit.aes=F,
                     shape=3, colour = 'gray25') +
          geom_point(x = cVWFA152[3], y = cVWFA152[2],inherit.aes=F,
                     shape=3, colour = 'gray25') +
          geom_point(x = pVWFA152[3], y = pVWFA152[2],inherit.aes=F,
                     shape=3, colour = 'gray25') +
          # geom_point(data = dfTractTh[dfTractTh$TRACT=='vof',],
          #             aes(x,y,color=probColor,fill=probColor),size=1.2) +
          # geom_point(data = dfTractTh[dfTractTh$TRACT=='parc',],
          #           aes(x,y,color=probColor,fill=probColor),size=1.2) +
          scale_colour_manual(name="", values = plotValues) +
          geom_point(aes(size=csize,colour=CLUSTER, fill=CLUSTER), alpha=0.2) +
          scale_size_identity() +
          geom_label(aes(x=czmean,y=cymean, label=TYPE),nudge_x = -.5, nudge_y = 1,
                     colour='black', size = 3, inherit.aes=F) +
          # geom_point(aes(size=cTmean),colour='black', fill='black', alpha = 0.7) +
          geom_point(aes(size=cTmean),colour='black', fill='black', alpha = 0.7) +
          xlab('MNI Z') + ylab('MNI Y') +
          scale_x_continuous(breaks=seq(-16,-4,2), limits= c(-17,-2)) +
          scale_y_continuous(breaks=seq(-80,-45,5), limits= c(-80,-44)) +
          # ggtitle(paste0('Clustered Contrasts, Tmin: ', Tmin)) +
          # facet_grid(TRT ~ DESIGN) +
          facet_grid(. ~ DESIGN) +
          coord_equal() +
          theme_bw() +
          theme(text  =element_text(size=8),
                legend.position="none",
                axis.text=element_text(size=8),
                # legend.position=c(.90, .90),
                axis.title=element_text(size=10,face="bold"),
                title=element_text(size=10, face="bold"),
                strip.text = element_text(size = 10, angle = 0,face="bold"),
                # panel.border = element_blank(),
                # panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
                # axis.line = element_line(colour = "black")
          )
        # DO NOT PLOT IT, SAVE IT AND RETURN IT
        # listToReturn[['clusterPlotZY']] = clusterPlotZY
        # # Plot: print on screen, high q vectorial, low q png
        # print(clusterPlotZY)
        # # svg(filename = file.path(PathFiguresFolder, paste0("clusterPlotZY_lq_",ver,".svg")),
        # #                 width = 18, height = 12, pointsize = 10)
        # # print(clusterPlot)
        # png(filename = file.path(PathFiguresFolder, paste0("clusterPlotZY_lq_",ver,".png")),
        #     width=400,height=550,units='px',pointsize=8,bg='white',res=90)
        # print(clusterPlotZY)
        # dev.off()
      }
    } # Cierra versionSinTractos
    versionConTractos = T
    if(versionConTractos){
      createPlotsXY = T
      if(createPlotsXY){
        clusterPlotXYcon = ggplot( data=df[df$DESIGN=='block',], aes(x=cymean,y=cxmean)) +
          geom_point(data = dfTractTh[dfTractTh$TRACT=='vof',],
                     aes(y,x,color=probColor,fill=probColor),size=1.2) +
          geom_point(data = dfTractTh[dfTractTh$TRACT=='parc',],
                     aes(y,x,color=probColor,fill=probColor),size=1.2) +
          # annotate('text', x = aVWFA152[1]-(lado/3), y = aVWFA152[2]+(lado/3),
          #          label='aVWFA',colour = 'gray25') +
          # annotate('text', x = cVWFA152[1]-(lado/3), y = cVWFA152[2]-(lado/3),
          #          label='cVWFA',colour = 'gray25') +
          # annotate('text', x = pVWFA152[1]-(lado/3), y = pVWFA152[2]-(lado/3),
          #          label='pVWFA',colour = 'gray25') +
          # geom_point(x = aVWFA152[1], y = aVWFA152[2],inherit.aes=F,
          #            shape=3, colour = 'gray25') +
          # geom_point(x = cVWFA152[1], y = cVWFA152[2],inherit.aes=F,
          #                shape=3, colour = 'gray25') +
          # geom_point(x = pVWFA152[1], y = pVWFA152[2],inherit.aes=F,
        #                shape=3, colour = 'gray25') +
        scale_colour_manual(name="", values = plotValues) +
          # geom_point(aes(size=csize, fill=CLUSTER,colour=CLUSTER), alpha=0.4) +
          geom_point(aes(size=csize), shape=1, colour='gray25', alpha=0.4) +
          scale_size_identity() +
          # geom_label(aes(x=cymean,y=cxmean, label=TYPE),nudge_x = 1, nudge_y = -4,
          #             colour='black', size = 3, inherit.aes=F) +
          # geom_label_repel(aes(x=cymean,y=cxmean, label=LABELTYPE),
          geom_label_repel(data=dft[dft$DESIGN=='block',], aes(x=yTYPE,y=xTYPE, label=LABELTYPE),
                           size = 3, inherit.aes=F, colour='black') +
          geom_point(aes(size=cTmean),colour='black', fill='black', alpha = 0.7) +
          xlab('MNI Y') + ylab('MNI X') +
          scale_y_reverse(breaks=seq(-60,0,10), limits= c(0,-60)) +
          scale_x_reverse(breaks=seq(-100,-10,10), limits= c(-10,-102)) +
          # ggtitle(paste0('Clustered clusters, Tmin: ', Tmin)) +
          # facet_grid(TRT ~ DESIGN) +
          # facet_grid(. ~ DESIGN) +
          coord_equal() +
          theme_bw() +
          theme(text  =element_text(size=12),
                legend.position="none",
                axis.text=element_text(size=12),
                # legend.position=c(.90, .90),
                axis.title=element_text(size=14,face="bold"),
                title=element_text(size=12, face="bold"),
                strip.text = element_text(size = 16, angle = 0,face="bold"),
                # panel.border = element_blank(),
                # panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
                # axis.line = element_line(colour = "black")
          )
        # DO NOT PLOT IT, SAVE IT AND RETURN IT
        listToReturn[['clusterPlotXYcon']] = clusterPlotXYcon
        # # Plot: print on screen, high q vectorial, low q png
        # # print(clusterPlotXY)
        # svg(filename = file.path(PathFigsFldrCap7,'tractContrastClusterAverages',
        #                          paste0("tractCLUSTERPlotXY_",ver,".svg")),
        #     width = 5, height = 3, pointsize = 12)
        # print(clusterPlotXY)
        # dev.off()
        # png(filename = file.path(PathFigsFldrCap7,'tractContrastClusterAverages',
        #                          paste0("tractCLUSTERPlotXY_",ver,".png")),
        #     width=5,height=3,units='in',pointsize=12,bg='white',res=300)
        # print(clusterPlotXY)
        # dev.off()
      }
      createPlotsZY = F
      if(createPlotsZY){
        clusterPlotZYcon = ggplot( data=df, aes(x=czmean,y=cymean)) +
          geom_point(data = dfTractTh[dfTractTh$TRACT=='vof',],
                     aes(z,y,color=probColor,fill=probColor),size=1.2) +
          geom_point(data = dfTractTh[dfTractTh$TRACT=='parc',],
                     aes(z,y,color=probColor,fill=probColor),size=1.2) +
          scale_colour_manual(name="", values = plotValues) +
          geom_point(aes(size=csize, fill=CLUSTER,colour=CLUSTER), alpha=0.4) +
          scale_size_identity() +
          geom_point(aes(x=czmean,y=cymean,size=cTmean),
                     colour='black',fill='black', alpha=0.6,inherit.aes=F) +
          annotate('text', x = aVWFA152[3]-(lado/3), y = aVWFA152[2]+(lado/3),
                   label='aVWFA',colour = 'gray25') +
          annotate('text', x = cVWFA152[3]-(lado/3), y = cVWFA152[2]-(lado/3),
                   label='cVWFA',colour = 'gray25') +
          annotate('text', x = pVWFA152[3]-(lado/3), y = pVWFA152[2]-(lado/3),
                   label='pVWFA',colour = 'gray25') +
          geom_point(x = aVWFA152[3], y = aVWFA152[2],inherit.aes=F,
                     shape=3, colour = 'gray25') +
          geom_point(x = cVWFA152[3], y = cVWFA152[2],inherit.aes=F,
                     shape=3, colour = 'gray25') +
          geom_point(x = pVWFA152[3], y = pVWFA152[2],inherit.aes=F,
                     shape=3, colour = 'gray25') +
          geom_label(aes(x=czmean,y=cymean, label=TYPE),nudge_x = -.5, nudge_y = 1,
                     colour='black', size = 3, inherit.aes=F) +
          geom_point(aes(size=cTmean),colour='black', fill='black', alpha = 0.7) +
          xlab('MNI Z') + ylab('MNI Y') +
          scale_x_continuous(breaks=seq(-30,30,10), limits= c(-30,30)) +
          scale_y_continuous(breaks=seq(-100,-10,10), limits= c(-102,-10)) +
          # ggtitle(paste0('Clustered clusters, Tmin: ', Tmin)) +
          # facet_grid(TRT ~ DESIGN) +
          facet_grid(. ~ DESIGN) +
          coord_equal() +
          theme_bw() +
          theme_bw() +
          theme(text  =element_text(size=12),
                legend.position="none",
                axis.text=element_text(size=12),
                # legend.position=c(.90, .90),
                axis.title=element_text(size=14,face="bold"),
                title=element_text(size=12, face="bold"),
                strip.text = element_text(size = 16, angle = 0,face="bold"),
                # panel.border = element_blank(),
                # panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
                # axis.line = element_line(colour = "black")
          )
        # DO NOT PLOT IT, SAVE IT AND RETURN IT
        listToReturn[['clusterPlotZYcon']] = clusterPlotZYcon
        # # Plot: print on screen, high q vectorial, low q png
        # # print(clusterPlotXY)
        # svg(filename = file.path(PathFigsFldrCap7,'tractContrastClusterAverages',
        #                          paste0("clusterPlotZY_",ver,".svg")),
        #     width = 6, height = 5, pointsize = 12)
        # print(clusterPlotZY)
        # dev.off()
      }
    } # Cierra versionConTractos

return(listToReturn)
}

# Perform all the anova analysis reported in the paper.
# Does not return anything, just show it in screen and store it in the R Notebook and version it with github
myAnova <- function(p) {

    # Read variables back
    dtmp                 = p[["dtmp"]]
    dependentVariable    = p[["dependentVariable"]]
    independentVariables = p[["independentVariables"]]

    # Choose between LME or LMER. Sometime one fails and the other one works.
    usarLME  = p[["usarLME"]]
    usarLMER = p[["usarLMER"]]

    # p.adjust.methods
    # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
    pAdjust = "fdr"




    hacerLeveneTest = F


    # **** EDIT THIS TO CHANGE ANALISIS (now the analysis is changed from the p parameters) ****
    grepId = '^SUBJECT$'
    # grepY   = '^fsy$'  # Select the ind. var. (las opciones son x, y, z, T, o fsx, fsy...)
    grepY   = paste0('^', dependentVariable, '$')

    if (length(independentVariables) == 1){
        grepXw1 = paste0('^', independentVariables[1], '$')
        grepXw2 = paste0('^', independentVariables[1], '$')
        grepXw3 = paste0('^', independentVariables[1], '$')
        # LME
        hacerConLME1Xw = T
        hacerConLME2Xw = F
        hacerConLME3Xw = F
        hacerPostHocAnalysisSinInteraccion1Xw = T
        hacerPostHocAnalysisSinInteraccion2Xw = F
        hacerPostHocAnalysis3Xw = F
        hacerPostHocAnalysisConInteraccion2Xw = F
        # LMER
        hacerConLMER1Xw = T
        hacerConLMER2Xw = F
        hacerConLMER3Xw = F
        LMERhacerPostHocAnalysisSinInteraccion1Xw = T
        LMERhacerPostHocAnalysisSinInteraccion2Xw = F
        LMERhacerPostHocAnalysisConInteraccion2Xw = F
        LMERhacerPostHocAnalysisSinInteraccion3Xw = F
    } else if (length(independentVariables) == 2){
        grepXw1 = paste0('^', independentVariables[1], '$')
        grepXw2 = paste0('^', independentVariables[2], '$')
        grepXw3 = paste0('^', independentVariables[2], '$')
        # LME
        hacerConLME1Xw = F
        hacerConLME2Xw = T
        hacerConLME3Xw = F
        hacerPostHocAnalysisSinInteraccion1Xw = F
        hacerPostHocAnalysisSinInteraccion2Xw = T
        hacerPostHocAnalysis3Xw = F
        hacerPostHocAnalysisConInteraccion2Xw = T
        # LMER
        hacerConLMER1Xw = F
        hacerConLMER2Xw = T
        hacerConLMER3Xw = F
        LMERhacerPostHocAnalysisSinInteraccion1Xw = F
        LMERhacerPostHocAnalysisSinInteraccion2Xw = T
        LMERhacerPostHocAnalysisConInteraccion2Xw = T
        LMERhacerPostHocAnalysisSinInteraccion3Xw = F
    } else if (length(independentVariables) == 3){
      grepXw1 = paste0('^', independentVariables[1], '$')
      grepXw2 = paste0('^', independentVariables[2], '$')
      grepXw3 = paste0('^', independentVariables[3], '$')

      # LME
      hacerConLME1Xw = F
      hacerConLME2Xw = F
      hacerConLME3Xw = T
      hacerPostHocAnalysisSinInteraccion1Xw = F
      hacerPostHocAnalysisSinInteraccion2Xw = F
      hacerPostHocAnalysis3Xw = T
      hacerPostHocAnalysisConInteraccion2Xw = F
      # LMER
      hacerConLMER1Xw = F
      hacerConLMER2Xw = F
      hacerConLMER3Xw = T
      LMERhacerPostHocAnalysisSinInteraccion1Xw = F
      LMERhacerPostHocAnalysisSinInteraccion2Xw = F
      LMERhacerPostHocAnalysisConInteraccion2Xw = F
      LMERhacerPostHocAnalysisSinInteraccion3Xw = T
    } else {
        cat(paste0('ERROR: ', length(independentVariables), ' is not an option for the number of independentVariables (should be 1, 2, 3)'))
    }
    # **** END EDIT THIS TO CHANGE ANALISIS ****


    crearDatos = T
    if(crearDatos) {
      # str(dtmp)
      # Do droplevels
      dtmp$SUBJECT = droplevels(dtmp$SUBJECT)
      dtmp$TRT = droplevels(dtmp$TRT)
      dtmp$DESIGN = droplevels(dtmp$DESIGN)
      dtmp$GROUP = droplevels(dtmp$GROUP)
      dtmp$ROI = droplevels(dtmp$ROI)
      # str(dtmp)

      # Asign names so that we do not need to edit the code below
      d1 = dtmp
      d1 = na.omit(d1)      ######## Check and think about this always
      names(d1)[grep(grepId , names(d1))] = 'id'
      names(d1)[grep(grepXw1, names(d1))] = 'Xw1'
      names(d1)[grep(grepXw2, names(d1))] = 'Xw2'
      names(d1)[grep(grepXw3, names(d1))] = 'Xw3'
      names(d1)[grep(grepY  , names(d1))] = 'Y'
      # str(d1)
    }



    # Hago el fit y compruebo si hay interaccion
    # usarLME = F
    if (usarLME){
      # hacerConLME1Xw = F
      if(hacerConLME1Xw){
        # Look at the averages first
        D1 = as.data.table(d1)
        D1[,media := mean(Y, na.rm =TRUE), by=list(Xw1)]
        D1[,sd :=sd(Y, na.rm =TRUE), by=list(Xw1)]
        print(unique(D1[, list(Xw1,media,sd)]))

        # Calculate the effect size of the difference
        # cohen.d(d1$Y,d1$Xw1)
        # cohen.d(d1$T,d1$Xw1)

        lmeFit = lme(Y ~ Xw1, random = ~ 1 | id, data=d1)
        print(anova(lmeFit))
      }
      # hacerConLME2Xw = F
      if(hacerConLME2Xw){
        # Look at the averages first
        D1 = as.data.table(d1)
        D1[,media := mean(Y, na.rm =TRUE), by=list(Xw1,Xw2)]
        D1[,sd :=sd(Y, na.rm =TRUE), by=list(Xw1,Xw2)]
        print(unique(D1[, list(Xw1,Xw2,media,sd)]))

        # Calculate the effect size of the difference
        # TEST
        # cohen.d(d1$Y[d1$Xw2=='TEST'],d1$Xw1[d1$Xw2=='TEST'])
        # cohen.d(d1$T[d1$Xw2=='TEST'],d1$Xw1[d1$Xw2=='TEST'])
        # cohen.d(d1$Y[d1$Xw2=='RETEST'],d1$Xw1[d1$Xw2=='RETEST'])
        # cohen.d(d1$T[d1$Xw2=='RETEST'],d1$Xw1[d1$Xw2=='RETEST'])


        # LME
        lmeReduced <- lme(Y ~ 1, random=list(id=pdBlocked(list(~1,
                                                               pdCompSymm(~Xw1-1), pdCompSymm(~Xw2-1)))),method="ML", data=d1)
        lmeFit = lme(Y ~ Xw1+Xw2, random=list(id=pdBlocked(list(~1,
                                                                pdCompSymm(~Xw1-1), pdCompSymm(~Xw2-1)))),method="ML", data=d1)
        lmeFitInteracc = lme(Y ~ Xw1*Xw2,random=list(id=pdBlocked(list(~1,
                                                                       pdCompSymm(~Xw1-1), pdCompSymm(~Xw2-1)))),method="ML", data=d1)
        print(anova(lmeFit))
        print(anova(lmeReduced, lmeFit));
        print(anova(lmeFitInteracc))
        print(anova(lmeReduced, lmeFitInteracc));
        print(anova(lmeFit, lmeFitInteracc));

      }
      # hacerConLME3Xw = T
      if(hacerConLME3Xw){
        # Look at the averages first
        D1 = as.data.table(d1)
        D1[,media := mean(Y, na.rm =TRUE), by=list(Xw1,Xw2,Xw3)]
        D1[,sd :=sd(Y, na.rm =TRUE), by=list(Xw1,Xw2,Xw3)]
        print(unique(D1[, list(Xw1,Xw2,Xw3,media,sd)]))


        lmeFit = lme(Y ~ Xw1+Xw2+Xw3, random=list(id=pdBlocked(list(~1,
                                                                    pdCompSymm(~Xw1-1), pdCompSymm(~Xw2-1),
                                                                    pdCompSymm(~Xw3-1)))),method="ML", data=d1)
        lmeFitInteracc = lme(Y ~ Xw1*Xw2*Xw3,random=list(id=pdBlocked(list(~1,
                                                                           pdCompSymm(~Xw1-1), pdCompSymm(~Xw2-1),
                                                                           pdCompSymm(~Xw3-1)))),method="ML", data=d1)
        print(anova(lmeFit, lmeFitInteracc)); print(anova(lmeFitInteracc))
      }
      # hacerPostHocAnalysisSinInteraccion1Xw = F
      if(hacerPostHocAnalysisSinInteraccion1Xw) {
        # Post-hoc analysis SIN INTERACCION
        print(summary(glht(lmeFit, linfct=mcp(Xw1="Tukey"),test=adjusted(pAdjust))))
        # confint(contrInteracc)  # Y si queremos saber los confidence intervals
      }
      # hacerPostHocAnalysisSinInteraccion2Xw = F
      if(hacerPostHocAnalysisSinInteraccion2Xw) {
        # Post-hoc analysis SIN INTERACCION
        print(summary(glht(lmeFit, linfct=mcp(Xw1="Tukey", Xw2="Tukey"),
                           test=adjusted(pAdjust))))
        # confint(contrInteracc)  # Y si queremos saber los confidence intervals
      }
      # hacerPostHocAnalysis3Xw = T
      if(hacerPostHocAnalysis3Xw) {
        # Post-hoc analysis SIN INTERACCION
        print(summary(glht(lmeFit, linfct=mcp(Xw1="Tukey",Xw2="Tukey",Xw3="Tukey"),
                           test=adjusted(pAdjust))))
        # Conf interval
        # confint(glht(lmeFit, linfct=mcp(Xw1="Tukey",Xw2="Tukey",Xw3="Tukey"),
        #                    test=adjusted(pAdjust)))
      }
      # hacerPostHocAnalysisConInteraccion2Xw = F
      if(hacerPostHocAnalysisConInteraccion2Xw) {
        # Post-hoc analysis CON INTERACCION
        tmp <- expand.grid(Xw1 = unique(d1$Xw1), Xw2 = unique(d1$Xw2))
        X <- model.matrix(~ Xw1 * Xw2, data = tmp)
        K = myObtener2factorK(as.data.frame(d1), 'Xw1', 'Xw2')
        contrInteracc = glht(lmeFitInteracc, linfct = as.matrix(K) %*% X)
        print(summary(contrInteracc))
        # confint(contrInteracc)  # Y si queremos saber los confidence intervals
      }
    }

    # CON LMER
    if (usarLMER){
      # OJO, para effect size:
      # https://stat.ethz.ch/pipermail/r-sig-mixed-models/2011q2/016134.html
      # http://stats.stackexchange.com/questions/95054/how-to-get-an-overall-p-value-and-effect-size-for-a-categorical-factor-in-a-mi
      r2.corr.mer <- function(m) {lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
      summary(lmfit)$r.squared}
      # r2.corr.mer(fitFyes)
      # o usar esto:
      # 1-var(residuals(fitFyes))/(var(model.response(model.frame(fitFyes))))
      # hacerConLMER1Xw = F
      if(hacerConLMER1Xw){
        # Look at the averages first
        D1 = as.data.table(d1)
        D1[,media := mean(Y, na.rm =TRUE), by=list(Xw1)]
        D1[,sd :=sd(Y, na.rm =TRUE), by=list(Xw1)]
        print(unique(D1[, list(Xw1,media,sd)]))

        # Calculate the effect size of the difference
        # cohen.d(d1$Y,d1$Xw1)
        # cohen.d(d1$T,d1$Xw1)

        # Ahora con lmer
        fitF <- lmer(Y~Xw1+(1|id),data=d1,REML=F)
        print(anova(fitF))
        summary(fitF)
        # Y para lmer ahora buscamos una p, comparamos el modelo con la base
        fitR <- lmer(Y ~ 1 + (1|id), data=d1, REML = F)
        print(KRmodcomp(fitF, fitR))
        AICc(fitF)
        print(aictab(cand.set=list(fitR, fitF), modnames=c("restricted", "full"),
                     sort=FALSE, second.ord=FALSE))
        # y ahora solo nos toca hacer las comparaciones multiples
        cat(paste0('Adjusted R-sq: ', r2.corr.mer(fitF)))
      }
      # hacerConLMER2Xw = T
      if(hacerConLMER2Xw){
        D1 = as.data.table(d1)
        D1[,media := mean(Y, na.rm =TRUE), by=list(Xw1,Xw2)]
        D1[,sd :=sd(Y, na.rm =TRUE), by=list(Xw1,Xw2)]
        print(unique(D1[, list(Xw1,Xw2,media,sd)]))


        # Ahora con lmer
        fitR <- lmer(Y ~ 1 + (1|id), data=d1, REML = F)
        fitFno <- lmer(Y~Xw1+Xw2+(1|id)+(1|Xw1:id)+(1|Xw2:id),data=d1,REML=F)
        fitFyes <- lmer(Y~Xw1*Xw2+(1|id)+(1|Xw1:id)+(1|Xw2:id),data=d1,REML=F)
        # print(anova(fitFno, fitFyes))
        #
        # # Tenemos que usar el modelo con la interaccion, a ver si puedo calcular las p-s
        # print(anova(fitFyes))
        # Y para lmer ahora buscamos una p, comparamos el modelo con la base
        print(KRmodcomp(fitR, fitFno))
        print(KRmodcomp(fitR, fitFyes))
        print(KRmodcomp(fitFno, fitFyes))
        # AICc(fitF)
        # print(aictab(cand.set=list(fitR, fitF), modnames=c("restricted", "full"),
        #        sort=FALSE, second.ord=FALSE))
        # y ahora solo nos toca hacer las comparaciones multiples
        cat(paste0('\nAdjusted R-sq para fitFno: ', r2.corr.mer(fitFno), '\n'))
        cat(paste0('\nAdjusted R-sq para fitFyes: ', r2.corr.mer(fitFyes), '\n'))
      }
      # hacerConLMER3Xw = F
      if(hacerConLMER3Xw){
        # Look at the averages first
        D1 = as.data.table(d1)
        D1[,media := mean(Y, na.rm =TRUE), by=list(Xw1,Xw2,Xw3)]
        D1[,sd :=sd(Y, na.rm =TRUE), by=list(Xw1,Xw2,Xw3)]
        print(unique(D1[, list(Xw1,Xw2,Xw3,media,sd)]))
        # Ahora con lmer
        fitR <- lmer(Y ~ 1 + (1|id), data=d1, REML = F)
        fitFno <- lmer(Y~Xw1+Xw2+Xw3+(1|id)+(1|Xw1:id)+(1|Xw2:id)+(1|Xw3:id),data=d1,REML=F)
        fitFyes <- lmer(Y~Xw1*Xw2*Xw3+(1|id)+(1|Xw1:id)+(1|Xw2:id)+(1|Xw3:id),data=d1,REML=F)
        # print(anova(fitFno, fitFyes))
        #
        # # Tenemos que usar el modelo con la interaccion, a ver si puedo calcular las p-s
        # print(anova(fitFyes))
        # Y para lmer ahora buscamos una p, comparamos el modelo con la base
        print(KRmodcomp(fitR, fitFno))
        print(KRmodcomp(fitR, fitFyes))
        print(KRmodcomp(fitFno, fitFyes))
        # AICc(fitF)
        # print(aictab(cand.set=list(fitR, fitF), modnames=c("restricted", "full"),
        #        sort=FALSE, second.ord=FALSE))
        # y ahora solo nos toca hacer las comparaciones multiples
        cat(paste0('\nAdjusted R-sq para fitFno: ', r2.corr.mer(fitFno), '\n'))
        cat(paste0('\nAdjusted R-sq para fitFyes: ', r2.corr.mer(fitFyes), '\n'))
      }
      # LMERhacerPostHocAnalysisSinInteraccion1Xw = F
      if(LMERhacerPostHocAnalysisSinInteraccion1Xw) {
        # Post-hoc analysis SIN INTERACCION
        print(summary(glht(fitF, linfct=mcp(Xw1="Tukey"),
                           test=adjusted(pAdjust))))
        # confint(contrInteracc)  # Y si queremos saber los confidence intervals
      }
      # LMERhacerPostHocAnalysisSinInteraccion2Xw = T
      if(LMERhacerPostHocAnalysisSinInteraccion2Xw) {
        # Post-hoc analysis SIN INTERACCION
        print(summary(glht(fitFno, linfct=mcp(Xw1="Tukey", Xw2="Tukey"),
                           test=adjusted(pAdjust))))
        # confint(contrInteracc)  # Y si queremos saber los confidence intervals
      }
      # LMERhacerPostHocAnalysisConInteraccion2Xw = T
      if(LMERhacerPostHocAnalysisConInteraccion2Xw) {
        # Post-hoc analysis CON INTERACCION
        tmp <- expand.grid(Xw1 = unique(d1$Xw1), Xw2 = unique(d1$Xw2))
        X <- model.matrix(~ Xw1 * Xw2, data = tmp)
        K = myObtener2factorK(as.data.frame(d1), 'Xw1', 'Xw2')
        contrInteracc = glht(fitFyes, linfct = as.matrix(K) %*% X)
        print(summary(contrInteracc))
        # confint(contrInteracc)  # Y si queremos saber los confidence intervals
      }
      # LMERhacerPostHocAnalysisSinInteraccion3Xw = F
      if(LMERhacerPostHocAnalysisSinInteraccion3Xw) {
        # Post-hoc analysis SIN INTERACCION
        print(summary(glht(fitFno, linfct=mcp(Xw1="Tukey",Xw2="Tukey",Xw3="Tukey"),
                           test=adjusted(pAdjust))))
      }
    }


    # hacerLeveneTest = F
    if(hacerLeveneTest){
      with(d1, leveneTest(Y, Xw2))
      with(d1, leveneTest(Y, interaction(Xw1, Xw2)))
      leveneTest(Y ~ Xw1*Xw2, data=d1)
      leveneTest(lm(Y ~ Xw1*Xw2, data=d1))
      leveneTest(Y ~ Xw1*Xw2, data=d1, center=mean)
      leveneTest(Y ~ Xw1*Xw2, data=d1, center=mean, trim=0.1)
    }


}


mySignalPercentChange <- function(trt, sujElim) {
  # trt = "TEST"
  # sujElim = sujElimArtRepair[['block']]

  resultsSPC = list()
  cat(trt, sujElim)

  # p.adjust.methods
  # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
  #   "fdr", "none")
  pAdjust = "fdr"



     # archivoSignalPercChange = '~/gDrive/BCBL/PROYECTOS/MINI/ANALYSIS/SignalPercentChange/block_GROUP60_block_OLD_sinMask_mpLOTS_roi-table.xlsx'
    # require(xlsx)
    # spcFILE = read.xlsx(archivoSignalPercChange, sheetIndex = 1)
    archivoSignalPercChange = '~/gDrive/BCBL/PROYECTOS/MINI/ANALYSIS/SignalPercentChange/TODO_block_GROUP60_block_OLD_sinMask_mpLOTS_roi-table_manualEdit.csv'
    spcFILE = read.csv(archivoSignalPercChange,header = TRUE,sep = ',')
    spcFILEs = spcFILE[, c('SUBJECT','CONTRAST','REGION','psc')]
    spcFILEs$REGION = mapvalues(spcFILEs$REGION, from = c("mLOTS", "pLOTS"), to = c("mOTS", "pOTS"))
    spcFILEs = spcFILEs[spcFILEs$SUBJECT %!in% sujElim, ]



    # require(data.table)
    spcdt = as.data.table(spcFILEs)

    # Remove outliers with SD = 2 (make them NA)
    multSD = 3
    spcdt[,out2sd := (psc < mean(psc)-multSD*sd(psc) | psc > mean(psc)+multSD*sd(psc)), by=.(CONTRAST, REGION)]
    spcdt$pscout2sd = spcdt$psc
    spcdt$pscout2sd[spcdt$out2sd == TRUE] = NA

    # Obtain the means
    spcdt[, MEAN := mean(pscout2sd, na.rm=TRUE), by=.(CONTRAST, REGION)]
    spcdt[, SD   :=   sd(pscout2sd, na.rm=TRUE), by=.(CONTRAST, REGION)]
    losSeis = c('RW-PW', 'RW-CS', 'RW-FF', 'RW-SD', 'RW-PS', 'RW-CB')
    # losSeis = c('PW-null', 'CS-null', 'FF-null', 'SD-null', 'PS-null', 'CB-null')
    spcdt = spcdt[spcdt$CONTRAST %in% losSeis,]
    spcdt$CONTRAST = factor(spcdt$CONTRAST, levels=losSeis)

    widedt =  data.table(CONTRAST=spcdt$CONTRAST[spcdt$REGION=='mOTS'],
                         mOTS=spcdt$pscout2sd[spcdt$REGION=='mOTS'],
                         pOTS=spcdt$pscout2sd[spcdt$REGION=='pOTS'])


    widedt[,ttpvalue := t.test(mOTS,pOTS,paired=TRUE)$p.value, by=CONTRAST]
    widedt[,ttestimate := t.test(mOTS,pOTS,paired=TRUE)$estimate, by=CONTRAST]
    uniquewidedt = unique(widedt[, .(CONTRAST, ttestimate,ttpvalue)])
    uniquewidedt[,SIGNIFICATIVO := ttpvalue < 0.05]
    widedt$CLUSTER[widedt$CONTRAST=='RW-PW'] = 'LEXICAL'
    widedt$CLUSTER[widedt$CONTRAST=='RW-CS'] = 'LEXICAL'
    widedt$CLUSTER[widedt$CONTRAST=='RW-FF'] = 'LEXICAL'
    widedt$CLUSTER[widedt$CONTRAST=='RW-PS'] = 'PERCEPTUAL'
    widedt$CLUSTER[widedt$CONTRAST=='RW-CB'] = 'PERCEPTUAL'
    widedt$CLUSTER[widedt$CONTRAST=='RW-SD'] = 'PERCEPTUAL'
    widedt[,gttpvalue := t.test(mOTS,pOTS,paired=TRUE)$p.value, by=CLUSTER]
    widedt[,gttestimate := t.test(mOTS,pOTS,paired=TRUE)$estimate, by=CLUSTER]
    guniquewidedt = unique(widedt[, .(CLUSTER, gttestimate,gttpvalue)])
    guniquewidedt[,SIGNIFICATIVO := gttpvalue < 0.05]










    spcdt$CLUSTER = NA
    spcdt$CLUSTER[spcdt$CONTRAST=='RW-PW'] = 'LEXICAL'
    spcdt$CLUSTER[spcdt$CONTRAST=='RW-CS'] = 'LEXICAL'
    spcdt$CLUSTER[spcdt$CONTRAST=='RW-FF'] = 'LEXICAL'
    spcdt$CLUSTER[spcdt$CONTRAST=='RW-PS'] = 'PERCEPTUAL'
    spcdt$CLUSTER[spcdt$CONTRAST=='RW-CB'] = 'PERCEPTUAL'
    spcdt$CLUSTER[spcdt$CONTRAST=='RW-SD'] = 'PERCEPTUAL'
    spcdt$CLUSTER[spcdt$CONTRAST=='PW-null'] = 'LEXICAL'
    spcdt$CLUSTER[spcdt$CONTRAST=='CS-null'] = 'LEXICAL'
    spcdt$CLUSTER[spcdt$CONTRAST=='FF-null'] = 'LEXICAL'
    spcdt$CLUSTER[spcdt$CONTRAST=='PS-null'] = 'PERCEPTUAL'
    spcdt$CLUSTER[spcdt$CONTRAST=='CB-null'] = 'PERCEPTUAL'
    spcdt$CLUSTER[spcdt$CONTRAST=='SD-null'] = 'PERCEPTUAL'
    spcdt$CLUSTER = factor(spcdt$CLUSTER, levels=c('LEXICAL','PERCEPTUAL'))





    # Ahora hacemos el anova
    # require(nlme)
    # require(pbkrtest)
    # library(multcomp)
    dtmp = copy(spcdt)
    grepId = '^SUBJECT$'
    grepXw3 = '^CONTRAST$'
    grepXw1 = '^REGION$'
    grepXw2 = '^CLUSTER$'
    grepY   = '^pscout2sd$'  # Select the ind. var.
    # **** END EDIT THIS TO CHANGE ANALISIS ****
    crearDatos = T
    if(crearDatos) {
      # Hacer asignaciones de nombres para luego no tener que tocar codigo
      d1 = dtmp
      # d1 = na.omit(d1)
      names(d1)[grep(grepId , names(d1))] = 'id'
      names(d1)[grep(grepXw1, names(d1))] = 'Xw1'
      names(d1)[grep(grepXw2, names(d1))] = 'Xw2'
      names(d1)[grep(grepXw3, names(d1))] = 'Xw3'
      names(d1)[grep(grepY  , names(d1))] = 'Y'
      d1$Xw1 = droplevels(d1$Xw1)
      d1$Xw2 = droplevels(d1$Xw2)
      d1$Xw3 = droplevels(d1$Xw3)
      str(d1)
    }
    # Hago el fit y compruebo si hay interaccion
    usarLME = F
    if (usarLME){
      hacerConLME1Xw = F
      if(hacerConLME1Xw){
        lmeFit = lme(Y ~ Xw1, random = ~ 1 | id, data=d1)
        print(anova(lmeFit))
      }
      hacerConLME2Xw = F
      if(hacerConLME2Xw){
        lmeReduced <- lme(Y ~ 1, random=list(id=pdBlocked(list(~1,
                                                               pdCompSymm(~Xw1-1), pdCompSymm(~Xw2-1)))),method="ML", data=d1)
        lmeFit = lme(Y ~ Xw1+Xw2, random=list(id=pdBlocked(list(~1,
                                                                pdCompSymm(~Xw1-1), pdCompSymm(~Xw2-1)))),method="ML", data=d1)
        lmeFitInteracc = lme(Y ~ Xw1*Xw2,random=list(id=pdBlocked(list(~1,
                                                                       pdCompSymm(~Xw1-1), pdCompSymm(~Xw2-1)))),method="ML", data=d1)
        print(anova(lmeFit))
        print(anova(lmeReduced, lmeFit));
        print(anova(lmeFitInteracc))
        print(anova(lmeReduced, lmeFitInteracc));
        print(anova(lmeFit, lmeFitInteracc));

      }
      hacerConLME3Xw = T
      if(hacerConLME3Xw){
        lmeFit = lme(Y ~ Xw1+Xw2+Xw3, random=list(id=pdBlocked(list(~1,
                                                                    pdCompSymm(~Xw1-1), pdCompSymm(~Xw2-1),
                                                                    pdCompSymm(~Xw3-1)))),method="ML", data=d1)
        lmeFitInteracc = lme(Y ~ Xw1*Xw2*Xw3,random=list(id=pdBlocked(list(~1,
                                                                           pdCompSymm(~Xw1-1), pdCompSymm(~Xw2-1),
                                                                           pdCompSymm(~Xw3-1)))),method="ML", data=d1)
        print(anova(lmeFit, lmeFitInteracc)); print(anova(lmeFitInteracc))
      }
      hacerPostHocAnalysisSinInteraccion1Xw = F
      if(hacerPostHocAnalysisSinInteraccion1Xw) {
        # Post-hoc analysis SIN INTERACCION
        print(summary(glht(lmeFit, linfct=mcp(Xw1="Tukey"),test=adjusted(pAdjust))))
        # confint(contrInteracc)  # Y si queremos saber los confidence intervals
      }
      hacerPostHocAnalysisSinInteraccion2Xw = F
      if(hacerPostHocAnalysisSinInteraccion2Xw) {
        # Post-hoc analysis SIN INTERACCION
        print(summary(glht(lmeFit, linfct=mcp(Xw1="Tukey", Xw2="Tukey"),
                           test=adjusted(pAdjust))))
        # confint(contrInteracc)  # Y si queremos saber los confidence intervals
      }
      hacerPostHocAnalysis3Xw = T
      if(hacerPostHocAnalysis3Xw) {
        # Post-hoc analysis SIN INTERACCION
        print(summary(glht(lmeFit, linfct=mcp(Xw1="Tukey",Xw2="Tukey",Xw3="Tukey"),
                           test=adjusted(pAdjust))))
        # Conf interval
        # confint(glht(lmeFit, linfct=mcp(Xw1="Tukey",Xw2="Tukey",Xw3="Tukey"),
        #                    test=adjusted(pAdjust)))
      }
      hacerPostHocAnalysisConInteraccion2Xw = F
      if(hacerPostHocAnalysisConInteraccion2Xw) {
        # Post-hoc analysis CON INTERACCION
        tmp <- expand.grid(Xw1 = unique(d1$Xw1), Xw2 = unique(d1$Xw2))
        X <- model.matrix(~ Xw1 * Xw2, data = tmp)
        K = myObtener2factorK(as.data.frame(d1), 'Xw1', 'Xw2')
        contrInteracc = glht(lmeFitInteracc, linfct = as.matrix(K) %*% X)
        print(summary(contrInteracc))
        # confint(contrInteracc)  # Y si queremos saber los confidence intervals
      }
    }
    usarLMER = T
    if (usarLMER){
      # CON LMER
      # OJO, para effect size:
      # https://stat.ethz.ch/pipermail/r-sig-mixed-models/2011q2/016134.html
      # http://stats.stackexchange.com/questions/95054/how-to-get-an-overall-p-value-and-effect-size-for-a-categorical-factor-in-a-mi
      r2.corr.mer <- function(m) {lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
      summary(lmfit)$r.squared}
      # r2.corr.mer(fitFyes)
      # o usar esto:
      # 1-var(residuals(fitFyes))/(var(model.response(model.frame(fitFyes))))
      hacerConLMER1Xw = F
      if(hacerConLMER1Xw){
        # Ahora con lmer
        fitF <- lmer(Y~Xw1+(1|id),data=d1,REML=F)
        print(anova(fitF))
        summary(fitF)
        # Y para lmer ahora buscamos una p, comparamos el modelo con la base
        fitR <- lmer(Y ~ 1 + (1|id), data=d1, REML = F)
        print(KRmodcomp(fitF, fitR))
        AICc(fitF)
        print(aictab(cand.set=list(fitR, fitF), modnames=c("restricted", "full"),
                     sort=FALSE, second.ord=FALSE))
        # y ahora solo nos toca hacer las comparaciones multiples
        cat(paste0('Adjusted R-sq: ', r2.corr.mer(fitF)))
      }
      hacerConLMER2Xw = T
      if(hacerConLMER2Xw){
        # Ahora con lmer
        fitR <- lmer(Y ~ 1 + (1|id), data=d1, REML = F)
        fitFno <- lmer(Y~Xw1+Xw2+(1|id)+(1|Xw1:id)+(1|Xw2:id),data=d1,REML=F)
        fitFyes <- lmer(Y~Xw1*Xw2+(1|id)+(1|Xw1:id)+(1|Xw2:id),data=d1,REML=F)
        # print(anova(fitFno, fitFyes))
        #
        # # Tenemos que usar el modelo con la interaccion, a ver si puedo calcular las p-s
        # print(anova(fitFyes))
        # Y para lmer ahora buscamos una p, comparamos el modelo con la base
        print(KRmodcomp(fitR, fitFno))
        print(KRmodcomp(fitR, fitFyes))
        print(KRmodcomp(fitFno, fitFyes))
        # AICc(fitF)
        # print(aictab(cand.set=list(fitR, fitF), modnames=c("restricted", "full"),
        #        sort=FALSE, second.ord=FALSE))
        # y ahora solo nos toca hacer las comparaciones multiples
        cat(paste0('\nAdjusted R-sq para fitFno: ', r2.corr.mer(fitFno), '\n'))
        cat(paste0('\nAdjusted R-sq para fitFyes: ', r2.corr.mer(fitFyes), '\n'))
      }
      hacerConLMER3Xw = F
      if(hacerConLMER3Xw){
        # Ahora con lmer
        fitR <- lmer(Y ~ 1 + (1|id), data=d1, REML = F)
        fitFno <- lmer(Y~Xw1+Xw2+Xw3+(1|id)+(1|Xw1:id)+(1|Xw2:id)+(1|Xw3:id),data=d1,REML=F)
        fitFyes <- lmer(Y~Xw1*Xw2*Xw3+(1|id)+(1|Xw1:id)+(1|Xw2:id)+(1|Xw3:id),data=d1,REML=F)
        # print(anova(fitFno, fitFyes))
        #
        # # Tenemos que usar el modelo con la interaccion, a ver si puedo calcular las p-s
        # print(anova(fitFyes))
        # Y para lmer ahora buscamos una p, comparamos el modelo con la base
        print(KRmodcomp(fitR, fitFno))
        print(KRmodcomp(fitR, fitFyes))
        print(KRmodcomp(fitFno, fitFyes))
        # AICc(fitF)
        # print(aictab(cand.set=list(fitR, fitF), modnames=c("restricted", "full"),
        #        sort=FALSE, second.ord=FALSE))
        # y ahora solo nos toca hacer las comparaciones multiples
        cat(paste0('\nAdjusted R-sq para fitFno: ', r2.corr.mer(fitFno), '\n'))
        cat(paste0('\nAdjusted R-sq para fitFyes: ', r2.corr.mer(fitFyes), '\n'))
      }
      LMERhacerPostHocAnalysisSinInteraccion1Xw = F
      if(LMERhacerPostHocAnalysisSinInteraccion1Xw) {
        # Post-hoc analysis SIN INTERACCION
        print(summary(glht(fitF, linfct=mcp(Xw1="Tukey"),
                           test=adjusted(pAdjust))))
        # confint(contrInteracc)  # Y si queremos saber los confidence intervals
      }
      LMERhacerPostHocAnalysisSinInteraccion2Xw = T
      if(LMERhacerPostHocAnalysisSinInteraccion2Xw) {
        # Post-hoc analysis SIN INTERACCION
        print(summary(glht(fitFno, linfct=mcp(Xw1="Tukey", Xw2="Tukey"),
                           test=adjusted(pAdjust))))
        # confint(contrInteracc)  # Y si queremos saber los confidence intervals
      }
      LMERhacerPostHocAnalysisConInteraccion2Xw = T
      if(LMERhacerPostHocAnalysisConInteraccion2Xw) {
        # Post-hoc analysis CON INTERACCION
        tmp <- expand.grid(Xw1 = unique(d1$Xw1), Xw2 = unique(d1$Xw2))
        X <- model.matrix(~ Xw1 * Xw2, data = tmp)
        K = myObtener2factorK(as.data.frame(d1), 'Xw1', 'Xw2')
        contrInteracc = glht(fitFyes, linfct = as.matrix(K) %*% X, alternative = c("greater"))
        print(summary(contrInteracc))
        # confint(contrInteracc)  # Y si queremos saber los confidence intervals
      }
      LMERhacerPostHocAnalysisSinInteraccion3Xw = F
      if(LMERhacerPostHocAnalysisSinInteraccion3Xw) {
        # Post-hoc analysis SIN INTERACCION
        print(summary(glht(fitFno, linfct=mcp(Xw1="Tukey",Xw2="Tukey",Xw3="Tukey"),
                           test=adjusted(pAdjust))))
      }
      hacerLeveneTest = F
      if(hacerLeveneTest){
        with(d1, leveneTest(Y, Xw2))
        with(d1, leveneTest(Y, interaction(Xw1, Xw2)))
        leveneTest(Y ~ Xw1*Xw2, data=d1)
        leveneTest(lm(Y ~ Xw1*Xw2, data=d1))
        leveneTest(Y ~ Xw1*Xw2, data=d1, center=mean)
        leveneTest(Y ~ Xw1*Xw2, data=d1, center=mean, trim=0.1)
      }
    }


    # MIRAMOS LAS MEDIAS, aunque haya dif significativas igual son ridiculas
    D1 = as.data.table(d1)

    # D1[,mediaT := mean(Y, na.rm =TRUE), by=list(Xw1)]
    # D1[,sdT :=sd(Y, na.rm =TRUE), by=list(Xw1)]
    # unique(D1[, list(Xw1,mediaT,sdT)])
    D1[,mediaT := mean(Y, na.rm=TRUE), by=list(Xw1,Xw2)]
    D1[,sdT :=sd(Y, na.rm=TRUE), by=list(Xw1,Xw2)]
    unique(D1[, list(Xw1,Xw2,mediaT,sdT)])
    # D1[,mediaT := mean(Y, na.rm =TRUE), by=list(Xw1,Xw2,Xw3)]
    # D1[,sdT :=sd(Y, na.rm =TRUE), by=list(Xw1,Xw2,Xw3)]
    # unique(D1[, list(Xw1,Xw2,Xw3,mediaT,sdT)])

    # model <- lm(mpg ~ hp + wt + hp:wt, data=mtcars)
    # library(effects)
    # plot(effect("hp:wt", model, list(wt=c(2.2,3.2,4.2))), multiline=TRUE)

    # END OF ANOVAS









    # spcmean = unique(spcdt[, .(CONTRAST, REGION,MEAN,SD)])
    spcmean = unique(D1[, list(Xw1,Xw2,mediaT,sdT)])
    names(spcmean) = c('REGION','GROUP','MEAN','SD')
    # N = 60  This was wrong
    N =length(unique(spcdt$SUBJECT))
    # cspcmean = spcmean[spcmean$CONTRAST %in% c('PW-null', 'CS-null', 'FF-null','SD-null', 'PS-null', 'CB-null'),]
    spcmean[, se := (SD/sqrt(N))]


    barplot <-
      ggplot(spcmean, aes(x=REGION, y=MEAN , colour=GROUP, fill=GROUP)) +
      geom_errorbar(aes(ymin=MEAN-se, ymax=MEAN+se), colour="black", width=.2) +
      geom_bar(position=position_dodge(), stat="identity",
               colour="black", # Use black outlines,
               size=.3) +      # Thinner lines
      # geom_point(position=pd, size=1, shape=21, fill="white") + # 21 is filled circle
      xlab("") +
      ylab("% signal change") +
      scale_fill_manual(name="", values = c("PERCEPTUAL"="#0072B2",
                                            "LEXICAL"="tomato3", # #D55E00",
                                            "TRANSITION"="#CC79A7")) +
      facet_grid(. ~ GROUP) +
      ggtitle("") +
      #expand_limits(y=0) +                        # Expand y range
      # scale_y_continuous(limits=c(0,4)) +         # Set tick every 4
      theme_bw() +
      theme(text  =element_text(size=10),
            legend.position="none",
            axis.text=element_text(size=10, face="bold"),
            axis.text.x = element_text(angle=0, hjust=0.5),
            # legend.position=c(.90, .90),
            axis.title=element_text(size=12,face="bold"),
            title=element_text(size=10, face="bold"),
            strip.text = element_text(size = 12, angle = 0,face="bold"),
            # panel.border = element_blank(),
            # panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            # axis.line = element_line(colour = "black")
      )

    resultsSPC[['barplot']] = barplot
    resultsSPC[['d']] = spcdt
    resultsSPC[['N']]       = length(unique(spcdt$SUBJECT))
    return(resultsSPC)


}



# LEXICAL DECISION FUNCTIONS
myReadLexicalDecision <- function() {
  # The lexical decision logs (the output from the Presentation tool) were first converted to csv with this python script.
  # ~/code/MINI/sReadLexicalDecisionLogs.py. This file does not do any conversion to the data.
  # All the data analysis is done in R in supportFunctions.R
  # I copied the extracted csv files to git to have always the original data for reproducibility purposes


    LDcsvdir = file.path(find_rstudio_root_file(), 'DATA','LD')
    DE = dir(LDcsvdir, pattern = '*log_ld.csv')
    LDGROUP = c(rep('ONCE', 35), rep('TWICE1', 30),
                rep('TWICE2', 30), 'TWICE1', 'TWICE2')
    LDPAIR  = c(rep('SOLO', 35),
                paste0('PAIR0', c(1:9)), paste0('PAIR', c(10:30)),
                paste0('PAIR0', c(1:9)), paste0('PAIR', c(10:30)),
                'PAIR31', 'PAIR31')
    LDTRT = c(rep('TEST', 64), rep('RETEST', 30), 'TEST', 'RETEST')
    LD = data.frame()
    # Read and data management

    for (di in seq_along(DE)){
      subject = substr(DE[di], 1, 4)
      # cat(subject)
      data = read.csv2(file.path(LDcsvdir,DE[di]), header = T, sep = ',')
      data = data[, c('Code', 'Type', 'RT')]
      names(data) = c('STIM', 'ACC', 'RT')
      con <- gsub('(.*)_(.*)_(.*)_(.*)', "\\2", data[,c('STIM')])
      con <- gsub('HighFRW', 'WH', con)
      con <- gsub('LowFRW' , 'WL', con)
      con <- gsub('CSCS'   , 'CS', con)
      con <- gsub('PWPW'   , 'PW', con)
      data$STIM = con
      data$STIM = as.factor(data$STIM)
      data$SUBJECT = as.factor(subject)
      data$GROUP = as.factor(LDGROUP[di])
      data$PAIR = as.factor(LDPAIR[di])
      data$TRT = as.factor(LDTRT[di])
      LD <- rbind(LD, data)
    }
    # Anadir al grupo al que pertenecen cada uno de los sujetos
    # He visto en el analisis que unos cuantos sujetos tienen un accuracy muy baja, < 10
    # Eso solo puede ser pq el contrabalanceo estuviera mal puesto
    # Apuntar sujetos y cambiar el orden
    changeOrderACC = c('S015','S059','S068','S074')
    hitInd       = LD$SUBJECT %in% changeOrderACC & LD$ACC == 'hit'
    incorrectInd = LD$SUBJECT %in% changeOrderACC & LD$ACC == 'incorrect'
    LD[hitInd, c('ACC')] = 'incorrect'
    LD[incorrectInd, c('ACC')] = 'hit'
    # ADemas estos sujetos tienen puntuacion muy baja, hacerlos NA
    #    SUBJECT    ACCSub  meanSub     sdSub
    # 2:    S029 72.909699 585.8545 198.96644
    # 4:    S067 61.872910 432.7926  92.19542
    # 7:    S097 61.204013 408.3044  87.43128
    ConvertToNA = c('S029','S067','S097')
    toNAInd = LD$SUBJECT %in% ConvertToNA
    LD[toNAInd, c('ACC')] = 'incorrect'
    # Convertir el reaction time de los incorrect a NA
    LD[LD$ACC=='incorrect',c('RT')] = NA


    # Ojo, hemos usado la info de abajo para hacer lo de arriba
    LDDT = as.data.table(LD)
    LDDT[ ,NSub        := .N, by=c('SUBJECT')]
    LDDT[ ,NSubStim    := .N, by=c('SUBJECT','STIM')]
    LDDT[ ,NStim       := .N, by=c('STIM')]
    LDDT[ ,ACCSub      := 100*sum(ACC=='hit', na.rm=TRUE)/NSub, by=c('SUBJECT')]
    LDDT[ ,ACCSubStim  := 100*sum(ACC=='hit', na.rm=TRUE)/NSubStim, by=c('SUBJECT', 'STIM')]
    LDDT[ ,ACCStim     := 100*sum(ACC=='hit', na.rm=TRUE)/NStim, by=c('STIM')]
    LDDT$RT = LDDT$RT / 10  #??Now it is in milliseconds
    inferiorLimit = 200
    LDDT$RT[LDDT$RT <= inferiorLimit] = NA
    # superiorLimit = 1500
    # LDDT$RT[LDDT$RT >= superiorLimit] = NA
    superiorLimitSD = 2
    LDDT$RT[LDDT$RT >= mean(LDDT$RT,na.rm=TRUE)+superiorLimitSD*sd(LDDT$RT,na.rm=TRUE)] = NA
    LDDT[ ,RTlog        := -log(1/RT)]
    LDDT[ ,meanSub      := mean(RT,na.rm=TRUE), by=c('SUBJECT')]
    LDDT[ ,sdSub        := sd(RT,na.rm=TRUE), by=c('SUBJECT')]
    LDDT[ ,meanStim     := mean(RT,na.rm=TRUE), by=c('STIM')]
    LDDT[ ,sdStim       := sd(RT,na.rm=TRUE), by=c('STIM')]
    LDDT[ ,meanSubStim  := mean(RT,na.rm=TRUE), by=c('SUBJECT','STIM')]
    LDDT[ ,sdSubStim    := sd(RT,na.rm=TRUE), by=c('SUBJECT','STIM')]

    # Todos aquellos ACC-s que sean 0, convertirlos a NA
    LDDT$ACCSub[LDDT$ACCSub==0] = NA
    LDDT$ACCSubStim[LDDT$ACCSubStim==0] = NA
    LDDT$ACCStim[LDDT$ACCStim==0] = NA


    LDDTSubStim = unique(LDDT[,list(SUBJECT,GROUP,STIM,ACCSubStim,meanSubStim)])
    # Algunos resultados son NaN-s, convertirlos a NAs.
    LDDTSubStim$meanSubStim[is.na(LDDTSubStim$meanSubStim)] = as.numeric(rep(NA,
                                                                             nrow(LDDTSubStim[is.na(LDDTSubStim$meanSubStim),.(meanSubStim)])))
    # Crear z-scores de cada uno de ellos y el combinado
    LDDTSubStim[ ,zRTSubStim   := scale(meanSubStim,TRUE,TRUE)]
    LDDTSubStim[ ,zACCSubStim  := scale(ACCSubStim ,TRUE,TRUE)]
    # OJO, lo he restado para el combinado, cuanto mas grande, mejor lector
    LDDTSubStim[ ,zCOMB        := zACCSubStim - zRTSubStim]
    return(LDDTSubStim)
}



myLDfMRIplots <- function(trt, sujElimArtRepair) {
  LDfMRIresults = list()

  # Read the information from the DATA folder
  dft = read.csv(file.path(find_rstudio_root_file(),
                            'DATA','LDfMRI',
                            paste0(trt,'_Data4BehavfMRIRegression.csv')),
                    header = T, sep = ',',dec = '.', fill = TRUE,
                    numerals = 'allow.loss', na.strings = 'NaN')

  # Check that we removed the subjects to be the same to the rest of the analyses
  #Remove the subjects identified as wrong
  dft = dft[dft$SUBJECT %!in% sujElimArtRepair, ]




  DT  = as.data.table(dft)
  DT[ ,N    := (.N - sum(is.na(LD))),by=c('LDcat','fMRIcat')]
  DT[ ,cor   :=my.cor.test(LD, spmT, 'estimate')  , by=c('LDcat','fMRIcat')]
  DT[ ,pval  :=my.cor.test(LD, spmT, 'p.value')  , by=c('LDcat','fMRIcat')]
  summaryResult = unique(DT[, list(LDcat,fMRIcat, N, cor, pval)])
  # Change order to annotate plots correctly
  summaryResult = summaryResult[c(2,1,4,3),]

  LDfMRIresults[['summaryResult']] = summaryResult





  # Now create a plot
  LDfMRIresults[['sctplot']] = ggplot(data=dft, aes(x=spmT, y=LD)) +
    geom_point(na.rm=TRUE, alpha = 0.4) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw() +
    coord_equal(ratio=2) +
    xlab('fMRI spmT values') +
    ylab('Lexical Decision RTs (zScore)') +
    annotate('text', x=8, y=1.15,label=paste0('r=',format(summaryResult$cor, digits=2)),hjust = 0) +
    annotate('text', x=8, y=0.8, label=paste0('p=',format(summaryResult$pval,digits=2)), hjust = 0) +
    annotate('text', x=8, y=0.45, label=paste0('N=',format(summaryResult$N,digits=2)), hjust = 0) +
    facet_grid(LDcat ~ fMRIcat) +
    scale_x_continuous(breaks = round(seq(min(dft$spmT), max(dft$spmT), by = 1),1)) +
    # scale_y_continuous(breaks = round(seq(min(dft$LD), max(dft$LD), by = 1),1)) +
    theme(text  =element_text(size=8),
          legend.position="none",
          axis.text=element_text(size=8,angle=0),
          # legend.position=c(.90, .90),
          axis.title=element_text(size=10,face="bold"),
          title=element_text(size=10, face="bold"),
          strip.text = element_text(size = 10, angle = 0,face="bold"),
          # panel.border = element_blank(),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
          # axis.line = element_line(colour = "black")
    )

  return(LDfMRIresults)





}




myqMRIplots <- function(trt, sujElimArtRepair) {
  qMRIresults = list()


  # FIX THIS, OBTAIN SUBJECT NAMES




  # Read the information from the DATA folder
  # mT1lex = read.csv(file.path(find_rstudio_root_file(), 'DATA','qMRI', paste0(trt, '_mT1lex4_T1.csv')),
  #                          header = F, sep = ',',dec = '.', fill = TRUE, numerals = 'allow.loss', na.strings = 'NaN')[,]
  # mT1per = read.csv(file.path(find_rstudio_root_file(), 'DATA','qMRI', paste0(trt, '_mT1per4_T1.csv')),
  #                          header = F, sep = ',',dec = '.', fill = TRUE, numerals = 'allow.loss', na.strings = 'NaN')[,]
  #
  # datos = data.frame(T1 = c(mT1lex, mT1per),
  #                    VWFA = factor(c(rep('mOTS',length(mT1lex)),rep('pOTS',length(mT1per))),
  #                                  levels = c('mOTS', 'pOTS')))

  datos = read.csv(file.path(find_rstudio_root_file(),
                             'DATA','qMRI', paste0(trt, '_mOTS_pOTS_T1.csv')),
                              header = T, sep = ',',dec = '.', fill = TRUE,
                              numerals = 'allow.loss', na.strings = 'NaN')
  datos = datos[datos$SUBJECT %!in% sujElimArtRepair, ]


  # In the first submission to pnas the outliers were removed, i don't know why, maybe I was testing something, I think they should not be removed
  # datos = myCleanOutliers(datos, 2, c(1))

  datos = na.omit(datos)
  qMRIresults[['ttest']] =t.test(datos$T1[datos$VWFA=='mOTS'],
                                 datos$T1[datos$VWFA=='pOTS'], paired = TRUE)
  qMRIresults[['effectSize']] = cohen.d(datos$T1, datos$VWFA, paired = T)



  qMRIresults[['violin']] = ggplot(data=datos, aes(x=VWFA, y=T1, color=VWFA)) +
    geom_violin(aes(fill = VWFA, alpha=0.8)) +
    geom_boxplot(width=.1,
                 outlier.shape=20,outlier.size=.75,outlier.stroke = 0.5, outlier.alpha=0.5) +
    # facet_grid(HEMI ~ Atlas) +
    theme_bw() +
    xlab('VWFA') +
    ylab('T1 relaxation [s]') +
    scale_fill_manual(name="", values = c("pOTS"="#0072B2", "mOTS"="tomato3")) +
    scale_color_manual(name="", values = c("pOTS"="#0072B2", "mOTS"="tomato3")) +
    theme(text  =element_text(size=8),
          legend.position="none",
          axis.text=element_text(size=8,angle=0),
          # legend.position=c(.90, .90),
          axis.title=element_text(size=10,face="bold"),
          title=element_text(size=10, face="bold"),
          strip.text = element_text(size = 10, angle = 0,face="bold"),
          # panel.border = element_blank(),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
          # axis.line = element_line(colour = "black")
    )

  datos_wide = data.frame(mT1lex = datos$T1[datos$VWFA=='mOTS'],
                          mT1per = datos$T1[datos$VWFA=='pOTS'])
  qMRIresults[['sctplot']] = ggplot(data=datos_wide, aes(x=mT1lex, y=mT1per)) +
    geom_point(na.rm=TRUE, alpha = 0.4) +
    stat_function(fun=function(x)x, geom='line',linetype="dashed") +
    theme_bw() +
    coord_equal() +
    xlab('mOTS T1 [s]') +
    ylab('pOTS T1 [s]') +
    annotate('text', x=1.5, y=1.25, label=paste0('N=',nrow(datos)/2)) +
    theme(text  =element_text(size=8),
          legend.position="none",
          axis.text=element_text(size=8,angle=0),
          # legend.position=c(.90, .90),
          axis.title=element_text(size=10,face="bold"),
          title=element_text(size=10, face="bold"),
          strip.text = element_text(size = 10, angle = 0,face="bold"),
          # panel.border = element_blank(),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
          # axis.line = element_line(colour = "black")
    )
  return(qMRIresults)
}



myAllspmTinVOTC <- function(trt, RWvsNull, indepContrasts, meanContrasts) {
  # trt = "TEST"
  # RWvsNull       = F
  # indepContrasts = F
  # meanContrasts  = T

  spmvotcResults = list()

    fMRIconts = c('VOT_block_RWvsCB','VOT_block_RWvsCS','VOT_block_RWvsFF',
                  'VOT_block_RWvsPS','VOT_block_RWvsPW','VOT_block_RWvsSD')
    TYPE = c('Perceptual','Lexical','Lexical','Perceptual','Lexical','Perceptual')
    spmTdir   = file.path(find_rstudio_root_file(), 'DATA','nogit_spmTvalues')
    fsSurfDir = file.path(find_rstudio_root_file(), 'DATA','nogit_fsaverage_surface')
    # Read files used for transformations
    lhwhite152 = read.csv(file.path(fsSurfDir,'lhwhite152.csv'), header = F,
                          sep = ',',dec = '.', fill = TRUE,
                          numerals = 'allow.loss', na.strings = 'NaN')

    if (RWvsNull){

      Nulltmp = read.csv(file.path(spmTdir, paste0(trt, '_VOT_block_RWvsNull.csv')),
                      header = F, sep = ',',dec = '.', fill = TRUE,
                      numerals = 'allow.loss', na.strings = 'NaN')
      NullMean = rowMeans(Nulltmp)
      spmDF = cbind(lhwhite152, NullMean)
      names(spmDF) = c('fsx','fsy','fsz','spmT')


      spmDFth = spmDF[spmDF$spmT >= 1.65,]
      spmDFth = spmDFth[spmDFth$fsy <= -30,]
      # nrow(spmDF)
      # nrow(spmDFth)
      # spmDFth[spmDFth$Contrast == 'RWvsNull',c('spmT')] = spmDFth[spmDFth$Contrast == 'RWvsNull',c('spmT')]/4
      # FFmax = spmDFth[spmDFth$spmT == max(spmDFth[spmDFth$Contrast == 'RWvsFF',c('spmT')]) &
      #                   spmDFth$Contrast == 'RWvsFF',
      #                 c('fsy', 'Contrast')]
      # Nullmax = spmDFth[spmDFth$spmT == max(spmDFth[spmDFth$Contrast == 'RWvsNull',c('spmT')]) &
      #                   spmDFth$Contrast == 'RWvsNull',
      #                   c('fsy','Contrast')]
      # PSmax = spmDFth[spmDFth$spmT == max(spmDFth[spmDFth$Contrast == 'RWvsPS',c('spmT')]) &
      #                   spmDFth$Contrast == 'RWvsPS',
      #                 c('fsy', 'Contrast')]
      Nullmax = spmDFth[spmDFth$spmT == max(spmDFth[,c('spmT')]),
                        c('fsy')]

      TNullPlot =  ggplot(data=spmDFth, aes(y=spmT,x=fsy,colour=spmT,fill=spmT)) +
        # TmapPlot =  ggplot(data=spmDFth, aes(y=spmT,x=fsx,colour=spmT,fill=spmT)) +
        # TmapPlot =  ggplot(data=spmDFth, aes(y=spmT,x=fsz,colour=spmT,fill=spmT)) +
        geom_point(size=2.5) +
        xlab('MNI Y') + ylab('SPM T value') +
        # xlab('MNI X') + ylab('SPM T value') +
        # xlab('MNI Z')   + ylab('SPM T value') +
        # ggtitle(paste0('Average RWvsNull spmTvalues >= 1.65 (p=0.05) in VOT')) +
        scale_colour_gradient2(low = 'yellow', mid = 'orange', high='red', midpoint = 4) +
        scale_fill_gradient2(low='yellow', mid = 'orange', high='red', midpoint = 4) +
        # facet_grid(. ~ Contrast ) +
        # coord_equal() +
        # coord_flip() +
        theme_bw() +
        scale_y_continuous(breaks=seq(2,10, 2), limits= c(1.6,11)) +
        scale_x_continuous(breaks=seq(-100,-30, 10), limits= c(-105,-30)) +
        # scale_y_continuous(breaks=seq(2,10, 2), limits= c(1.6,11)) +
        # scale_x_continuous(breaks=seq(-50,-10, 10), limits= c(-55,-10)) +
        # scale_y_continuous(breaks=seq(2,10, 2), limits= c(1.6,11)) +
        # scale_x_continuous(breaks=seq(-20,20, 5), limits= c(-20,20)) +
        theme(text  =element_text(size=8),
              legend.position="none",
              axis.text=element_text(size=8),
              # legend.position=c(.90, .90),
              axis.title=element_text(size=10,face="bold"),
              title=element_text(size=10, face="bold"),
              strip.text = element_text(size = 10, angle = 0,
                                        face="bold"))

      # print(TmapPlot)
      # svg(filename = file.path(svgbaseMINIimgDir, paste0("RWvsNull_spmTmap_Z.svg")),
      #                 width = 18, height = 12, pointsize = 10)
      # print(TmapPlot)
      # png(filename = file.path(pngbaseMINIimgDir, paste0("RWvsNull_spmTmap_Y.png")),
      #                 width = 800, height = 400, units = 'px', pointsize = 12, bg ='white',res=300)
      # print(TmapPlot)
      # dev.off()
      spmvotcResults[['TNullPlot']] = TNullPlot
    }



    # Create the data tables for the different contrast.
    # They will be used to plot them directly or to create the LEX and PER averages
    maximosY = list()
    maximosT = list()
    minimosT = list()
    A = data.frame()
    # TH  = 0  Not necessary, done in Matlab now
    for (i in 1:length(TYPE)){
      spmT = read.csv(file.path(spmTdir, paste0(trt, '_',fMRIconts[i],'.csv')),
                      header = FALSE, sep=',', dec = '.', fill = TRUE,
                      numerals = 'allow.loss', na.strings = 'NaN',
                      row.names = NULL)
      # spmT[spmT <= TH] = NA
      spmT = rowMeans(spmT, na.rm =TRUE)
      temp = lhwhite152
      names(temp) = c('fsx','fsy','fsz')
      temp$spmT = spmT
      spmDF = temp
      spmDF$Contrast = fMRIconts[i]
      spmDF$TYPE = TYPE[i]

      # Delete, not used, I think it was a test
      # spmDFth = spmDF[!is.na(spmDF$fsy),
      # spmDFth = spmDFth[spmDFth$spmT >= 1.65,]
      # spmDFth = spmDFth[spmDFth$fsy <= -30,]

      # Y ahora anadirlo al data frame global
      A = rbind(A, spmDF)
    }
    A$Contrast = as.factor(A$Contrast)
    A$TYPE = as.factor(A$TYPE)
    A = A[,c('TYPE', 'Contrast', 'fsy','spmT')]
    A$fsy = round(A$fsy,1)
    # str(A)
    B = as.data.table(A)
    # str(B)
    absmax <- function(x) { x[which.max( abs(x) )]}
    B[,spmTmaxY  := absmax(spmT),by=.(Contrast,fsy)]
    B[,spmTmeanY := mean(spmT, na.rm=TRUE), by=.(Contrast,fsy)]

    Cmax = unique(B[,.(TYPE, Contrast,fsy,spmTmaxY)])
    C = Cmax

    C = C[C$fsy        <= -40,]
    C = C[C$fsy        >= -80,]

    C$Contrast = factor(C$Contrast, levels=c("VOT_block_RWvsPS", "VOT_block_RWvsCB","VOT_block_RWvsSD",
                                             "VOT_block_RWvsFF","VOT_block_RWvsCS","VOT_block_RWvsPW"))

    if (indepContrasts){

      # TmapPlot =  ggplot(data=A, aes(y=spmT,x=fsy,colour=spmT,fill=spmT)) +
      TcontrastsPlot =  ggplot(data=C, aes(y=spmTmaxY,x=fsy,colour=Contrast,fill=Contrast)) +
        # TmapPlot =  ggplot(data=C, aes(y=spmTmeanY,x=fsy,colour=Contrast,fill=Contrast)) +
        geom_line(size = .8, alpha=.75 )+ #??, position='identity') +  #
        # geom_area(size = .8, alpha=0.4, position='identity') +  #
        # stat_smooth(method="lm", se=FALSE, fill=NA,formula=y ~ poly(x, 40, raw=TRUE)) +
        # geom_bin2d(bins = 100,alpha = 0.5) +  #
        xlab('MNI Y') + ylab('Average SPM T value per Vertex') +
        # scale_colour_gradient2(low = 'yellow', mid = 'orange', high='red', midpoint = 4) +
        # scale_fill_gradient2(low='yellow', mid = 'orange', high='red', midpoint = 4) +
        # facet_grid(TYPE ~ . , scales = 'free') +
        # facet_grid(. ~ TYPE ) +
        # coord_equal() +
        theme_bw() +
        scale_x_reverse() +
        # scale_y_continuous(breaks=seq(-2.5,2.5, .5), limits= c(-3,3)) +
        # scale_x_continuous(breaks=seq(-100,-30, 10), limits= c(-105,-30)) +
        # ylim(c(0.5,2.5)) +
        theme(text  =element_text(size=12),
              # legend.position="none",
              axis.text=element_text(size=12),
              legend.position=c(.90, .90),
              axis.title=element_text(size=14,face="bold"),
              title=element_text(size=14, face="bold"),
              strip.text = element_text(size = 14, angle = 0, face="bold"))
      spmvotcResults[['TcontrastsPlot']] = TcontrastsPlot

    }

    # Now do the same but average it over the LEXICAL and the PERCEPTUAL contrasts
    # We start from C
    if (meanContrasts){
      D = copy(C)
      D[, spmTmaxY2 := mean(spmTmaxY, na.rm=TRUE), by=.(TYPE, fsy)]
      D = unique(D[,.(TYPE, fsy, spmTmaxY2)])

      # Plot wothout binning
      # TmapPlot =  ggplot(data=D, aes(y=spmTmaxY2,x=fsy,colour=TYPE,fill=TYPE)) +
      #   geom_line(size = .8, alpha=.75 ) +
      #   xlab('MNI Y') + ylab('Average SPM T value per Vertex') +
      #   theme_bw() +
      #   scale_x_reverse() +
      #   theme(text  =element_text(size=12),
      #         axis.text=element_text(size=12),
      #         legend.position=c(.10, .90),
      #         axis.title=element_text(size=14,face="bold"),
      #         title=element_text(size=14, face="bold"),
      #         strip.text = element_text(size = 14, angle = 0, face="bold"))


      # Now do the same but try to bin the y values
      dfD    = as.data.frame(D)
      dfDper = as.data.frame(D[D$TYPE=='Perceptual',])
      dfDwl  = as.data.frame(D[D$TYPE=='Lexical',])
      dfDperbin = ddply(dfDper, .(cut( dfDper$fsy, ((80-40)*5))), summarize,
                        meanSPM = round(mean(spmTmaxY2), 2),
                        meanfsy = round(mean(fsy), 2))
      dfDwlbin = ddply(dfDwl, .(cut( dfDwl$fsy, ((80-40)*5))), summarize,
                       meanSPM = round(mean(spmTmaxY2), 2),
                       meanfsy = round(mean(fsy), 2))
      E = cbind(c(rep('PERCEPTUAL', nrow(dfDperbin)), rep('LEXICAL', nrow(dfDwlbin))),
                rbind(dfDperbin[,c('meanSPM', 'meanfsy')], dfDwlbin[,c('meanSPM', 'meanfsy')]))
      names(E)[1] = 'ContrastType'
      E$ContrastType = factor(E$ContrastType, levels = c('PERCEPTUAL', 'LEXICAL'))

      TmeanPlot =  ggplot(data=E, aes(y=meanSPM,x=meanfsy,colour=ContrastType,fill=ContrastType)) +
        geom_line(size = 1.2, alpha=.75 ) +
        xlab('MNI Y') + ylab('Average SPM T value per Vertex') +
        theme_bw() +
        coord_equal(ratio=20) +
        scale_x_reverse() +
        theme(text  =element_text(size=12),
              axis.text=element_text(size=12),
              legend.position=c(.15, .80),
              axis.title=element_text(size=14,face="bold"),
              title=element_text(size=14, face="bold"),
              strip.text = element_text(size = 14, angle = 0, face="bold"))

      spmvotcResults[['TmeanPlot']] = TmeanPlot
  }

return(spmvotcResults)
}


