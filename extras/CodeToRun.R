library(BiasInAcetaminophenStudies)

options(fftempdir = "s:/FFtemp")
connectionDetails <- createConnectionDetails(dbms = "pdw",
                                             server = Sys.getenv("PDW_SERVER"),
                                             port = Sys.getenv("PDW_PORT"))
cdmDatabaseSchema <- "cdm_cprd_v1017.dbo"
cohortDatabaseSchema <- "scratch.dbo"
cohortTable <- "mschuemi_epi_688_cohorts"
oracleTempSchema <- NULL
outputFolder <- "s:/BiasInAcetaminophenStudies"
maxCores <- parallel::detectCores()

execute(connectionDetails = connectionDetails,
        cdmDatabaseSchema = cdmDatabaseSchema,
        cohortDatabaseSchema = cohortDatabaseSchema,
        cohortTable = cohortTable,
        oracleTempSchema = oracleTempSchema,
        outputFolder = outputFolder,
        createCohorts = FALSE,
        runCohortMethod = FALSE,
        runCaseControl = FALSE,
        createPlotsAndTables = TRUE,
        generateReport = TRUE,
        maxCores = maxCores,
        blind = FALSE) 
