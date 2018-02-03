context("Survival analysis")

clinical <- read.table(text = "2549   NA ii  female
                                840   NA i   female
                                 NA 1204 iv    male
                                 NA  383 iv  female
                               1293   NA iii   male
                                 NA 1355 ii    male")
names(clinical) <- c("patient.days_to_last_followup", "patient.days_to_death",
                     "patient.stage_event.pathologic_stage", "patient.gender")

test_that("Test processing survival terms", {
    timeStart  <- "days_to_death"
    event      <- "days_to_death"
    formulaStr <- "patient.stage_event.pathologic_stage + patient.gender"
    survTerms <- processSurvTerms(clinical, censoring="right", event, timeStart,
                                  formulaStr=formulaStr)
    expect_is(survTerms, "list")
    expect_length(survTerms, 3)
    expect_is(survTerms$form, "formula")
    expect_equal(as.character(survTerms$form)[[3]], formulaStr)
    expect_is(survTerms$survTime, "data.frame")
    expect_true( all(survTerms$survTime[5:8] == clinical, na.rm=TRUE) )
    expect_is(survTerms$scale, "character")
    expect_equal(survTerms$scale, "days")
    
    # Check events
    event <- !is.na(clinical$patient.days_to_death)
    expect_equal(unique(survTerms$survTime$event[event]), 1)
    expect_equal(unique(survTerms$survTime$event[!event]), 0)
    
    # Check time
    time <- clinical$patient.days_to_death
    time[!event] <- clinical$patient.days_to_last_followup[!event]
    expect_equal(survTerms$survTime$time, time)
})

test_that("Test processing survival terms to fit a Cox PH model", {
    timeStart  <- "days_to_death"
    event      <- "days_to_death"
    formulaStr <- "patient.stage_event.pathologic_stage + patient.gender"
    expect_warning(
        survTerms <- processSurvTerms(clinical, censoring="right", event, 
                                      timeStart, formulaStr=formulaStr, 
                                      coxph = TRUE))
    expect_is(survTerms, "coxph")
    expect_equal(survTerms$n, nrow(clinical))
    
    event <- !is.na(clinical$patient.days_to_death)
    expect_equal(survTerms$nevent, sum(event))
})

test_that("Test survival difference between groups", {
    timeStart  <- "days_to_death"
    event      <- "days_to_death"
    formulaStr <- "patient.stage_event.pathologic_stage + patient.gender"
    survTerms <- processSurvTerms(clinical, censoring="right", event, timeStart,
                                  formulaStr=formulaStr)
    pvalue <- testSurvival(survTerms)
    expect_is(pvalue, "numeric")
    expect_equal(pvalue, 0.196, tolerance=1e-3)
})

test_that("Plot survival curves", {
    timeStart  <- "days_to_death"
    event      <- "days_to_death"
    formulaStr <- "patient.stage_event.pathologic_stage + patient.gender"
    survTerms <- processSurvTerms(clinical, censoring="right", event, timeStart,
                                  formulaStr=formulaStr)
    surv <- survfit(survTerms)
    pvalue <- testSurvival(survTerms)

    plot <- plotSurvivalCurves(surv, pvalue=pvalue)
    expect_is(plot, "highchart")
    expect_equal(plot$x$type, "chart")
    expect_equivalent(plot$x$hc_opts$yAxis[c("min", "max")], c(0, 1))
    expect_match(plot$x$hc_opts$subtitle$text, as.character(pvalue))
    expect_equal(plot$x$hc_opts$chart$zoomType, "xy")
    expect_length(plot$x$hc_opts$series, nrow(clinical))
})

test_that("Quantify optimal PSI cutoff", {
    timeStart  <- "days_to_death"
    event      <- "days_to_death"
    data       <- c(0.1, 0.2, 0.9, 1, 0.2, 0.6)
    
    opt <- optimalSurvivalCutoff(clinical, data, censoring="right", event, 
                                 timeStart)
    expect_is(opt, "list")
    expect_equal(cutoff <- opt$par, 0.618034, tol=1e-6)
    expect_equal(pvalue <- opt$value, 0.0269, tol=1e-4)
    expect_equal(opt$convergence, 0)
})

test_that("Plot survival curves separated by PSI cutoff", {
    cutoff <- 0.5
    data   <- c(0.1, 0.2, 0.9, 1, 0.2, 0.6)
    group  <- labelBasedOnCutoff(data, cutoff, "Inclusion levels")
    
    timeStart  <- "days_to_death"
    event      <- "days_to_death"
    survTerms <- processSurvTerms(clinical, censoring="right", event, timeStart,
                                  group=group)
    surv <- survfit(survTerms)
    pvalue <- testSurvival(survTerms)
    plot <- plotSurvivalCurves(surv, pvalue=pvalue)
    
    expect_is(plot, "highchart")
    expect_equal(plot$x$type, "chart")
    expect_equivalent(plot$x$hc_opts$yAxis[c("min", "max")], c(0, 1))
    expect_match(plot$x$hc_opts$subtitle$text, as.character(pvalue))
    expect_equal(plot$x$hc_opts$chart$zoomType, "xy")
    expect_length(plot$x$hc_opts$series, 2)
})

test_that("Fit a Cox PH model for PSI cutoff separation", {
    cutoff <- 0.5
    data   <- c(0.1, 0.2, 0.9, 1, 0.2, 0.6)
    group <- labelBasedOnCutoff(data, cutoff, "Inclusion levels")
    
    timeStart  <- "days_to_death"
    event      <- "days_to_death"
    expect_warning(
        survTerms <- processSurvTerms(clinical, censoring="right", event,
                                      timeStart, group=group, coxph=TRUE))
    
    expect_is(survTerms, "coxph")
    expect_equal(survTerms$n, nrow(clinical))
    
    event <- !is.na(clinical$patient.days_to_death)
    expect_equal(survTerms$nevent, sum(event))
})

test_that("Plot survival curves with no separation", {
    cutoff <- 0.6
    data   <- c(0.1, 0.2, 0.9, 1, 0.2, 0.6)/2
    group  <- labelBasedOnCutoff(data, cutoff, "Inclusion levels")
    
    timeStart  <- "days_to_death"
    event      <- "days_to_death"
    survTerms <- processSurvTerms(clinical, censoring="right", event, timeStart,
                                  group=group)
    surv <- survfit(survTerms)
    pvalue <- testSurvival(survTerms)
    plot <- plotSurvivalCurves(surv, pvalue=pvalue)
    
    expect_is(plot, "highchart")
    expect_equal(plot$x$type, "chart")
    expect_equivalent(plot$x$hc_opts$yAxis[c("min", "max")], c(0, 1))
    expect_match(plot$x$hc_opts$subtitle$text, "NA")
    expect_equal(plot$x$hc_opts$chart$zoomType, "xy")
    expect_length(plot$x$hc_opts$series, 1)
})