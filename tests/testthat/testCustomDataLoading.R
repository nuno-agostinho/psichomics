context("Functions to load and process user-provided data")

test_that("Prepare SRA metadata", {
    # SRA file with metadata
    sraMetadata <- NULL
    sraMetadata$header <- paste0(
        'Run,Assay Type,AvgSpotLen,BioProject,BioSample,Cell_Line,Cell_type,',
        'Center Name,Consent,DATASTORE filetype,DATASTORE provider,',
        'DATASTORE region,Experiment,genotype/variation,GEO_Accession,',
        'Instrument,LibraryLayout,LibrarySelection,LibrarySource,MBases,',
        'MBytes,Organism,passage,Platform,ReleaseDate,sample_acc,Sample Name,',
        'source_name,SRA Study')
    sraMetadata$line1 <- paste0(
        'SRR6368612,RNA-Seq,300,PRJNA422102,SAMN08164377,HT29,human colorectal',
        'cancer cells,GEO,public,sra,"gs,ncbi,s3","gs.US,ncbi.public,',
        's3.us-east-1",SRX3464014,wild type; control,GSM2885018,Illumina HiSeq',
        '2500,PAIRED,cDNA,TRANSCRIPTOMIC,4319,1676,Homo sapiens,16-19,',
        'ILLUMINA,2018-03-12T00:00:00Z,SRS2752018,GSM2885018,HT29_Control',
        'siRNA,SRP126561')
    sraMetadata$line2 <- paste0(
        'SRR6368613,RNA-Seq,300,PRJNA422102,SAMN08164376,HT29,human colorectal',
        'cancer cells,GEO,public,sra,"gs,ncbi,s3","gs.US,ncbi.public,',
        's3.us-east-1",SRX3464015,wild type; control,GSM2885019,Illumina HiSeq',
        '2500,PAIRED,cDNA,TRANSCRIPTOMIC,5451,2120,Homo sapiens,16-19,',
        'ILLUMINA,2018-03-12T00:00:00Z,SRS2752023,GSM2885019,HT29_Control',
        'siRNA,SRP126561')
    sraMetadata$line3 <- paste0(
        'SRR6368614,RNA-Seq,300,PRJNA422102,SAMN08164375,HT29,human colorectal',
        'cancer cells,GEO,public,sra,"gs,ncbi,s3","gs.US,ncbi.public,',
        's3.us-east-1",SRX3464016,wild type; control,GSM2885020,Illumina HiSeq',
        '2500,PAIRED,cDNA,TRANSCRIPTOMIC,5676,2345,Homo sapiens,16-19,',
        'ILLUMINA,2018-03-12T00:00:00Z,SRS2752019,GSM2885020,HT29_Control ',
        'siRNA,SRP126561')
    sraMetadata$line4 <- paste0(
        'SRR6368615,RNA-Seq,300,PRJNA422102,SAMN08164374,HT29,human colorectal',
        ' cancer cells,GEO,public,sra,"gs,ncbi,s3","gs.US,ncbi.public,',
        's3.us-east-1",SRX3464017,RNF6 knockdown,GSM2885021,Illumina HiSeq ',
        '2500,PAIRED,cDNA,TRANSCRIPTOMIC,6911,2704,Homo sapiens,16-19,',
        'ILLUMINA,2018-03-12T00:00:00Z,SRS2752020,GSM2885021,HT29_RNF6 siRNA,',
        'SRP126561')
    sraMetadata <- paste(unlist(sraMetadata), collapse="\n")
    
    info <- prepareSRAmetadata(sraMetadata, output=NULL)
    expect_identical(info$`Sample ID`, paste0("SRR636861", 2:5))
    expect_identical(unique(info$`Assay Type`), "RNA-Seq")
    expect_identical(unique(info$Cell_Line), "HT29")
    expect_identical(unique(info$`Center Name`), "GEO")
    expect_identical(unique(info$Consent), "public")
    expect_identical(unique(info$LibraryLayout), "PAIRED")
    expect_identical(unique(info$LibrarySelection), "cDNA")
    expect_identical(unique(info$LibrarySource), "TRANSCRIPTOMIC")
    expect_identical(unique(info$Platform), "ILLUMINA")
    expect_identical(unique(info$passage), "16-19")
    expect_identical(unique(info$Organism), "Homo sapiens")
})
