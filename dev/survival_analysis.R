# Load BRCA data from TCGA
data               <- loadTCGAdata(cohort="BRCA")[[1]]

brca               <- list()
brca$clinical      <- data$`Clinical data`
brca$sampleInfo    <- data$`Sample metadata`
brca$junctionQuant <- data$`Junction quantification (Illumina HiSeq)`

# Quantify alternative splicing
annot <- loadAnnotation(listSplicingAnnotations(assembly="hg19")[[1]])
psi   <- quantifySplicing(annot, junctionQuant)

# AS events in question
events <- c("SE_7_-_23571408_23562051_23561751_23561459_TRA2A",
            "SE_1_+_46072263_46072993_46074009_46078841_NASP",
            "SE_17_+_62510476_62512841_62512946_62515439_CEP95",
            "SE_16_-_2318056_2314761_2314574_2314332_RNPS1",
            "SE_1_-_228296850_228296722_228296656_228296019_MRPL55")

match <- getSubjectFromSample(colnames(psi), clinical, sampleInfo=sampleInfo)
eventData <- lapply(events, function(k) assignValuePerPatient(psi[k, ], match))
names(eventData) <- events

event <- timeStart <- "days_to_death"
optm  <- optimalSurvivalCutoff(clinical, eventData[[1]], censoring="right",
                               timeStart=timeStart, event=event)
testSurvivalCutoff(0.305, eventData[[2]], clinical=clinical,
                   timeStart=timeStart, event=event, censoring=censoring)

# Plot survival
plotSurvivalPvaluesByCutoff(
    clinical, eventData[[1]],
    timeStart=timeStart, event=event, censoring=censoring)
