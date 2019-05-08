# Install/Load required packages

if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
require("httr")

if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
require("jsonlite")

if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
require("data.table")

options(repos = BiocManager::repositories())
getOption("repos")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("VariantAnnotation", quietly = TRUE)) BiocManager::install("VariantAnnotation")
require("VariantAnnotation")
if (!requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
require("BSgenome.Hsapiens.1000genomes.hs37d5")

if (!file.exists("data/Txdb.sqlite")) {
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) BiocManager::install("VariantAnnotation")
  require("GenomicFeatures")
  #Txdb <- makeTxDbFromEnsembl(organism="Homo sapiens", release=75, circ_seqs=DEFAULT_CIRC_SEQS, server="ensembldb.ensembl.org", username="anonymous", password=NULL, port=0L, tx_attrib=NULL)
  #### Fix the 5' CDS Incomplete issue and regenerate the Txdb from scratch
  gtf <- fread("./data/Homo_sapiens.GRCh37.75.gtf", header = F, skip = "#", stringsAsFactors = F, na.strings = "")
  gtf2 <- gtf[V1 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT")][, ':='(geneid = gsub(".*(ENSG\\d+).*", "\\1", V9), transcriptid = gsub(".*(ENST\\d+).*", "\\1", V9), exonnumber = ifelse(V3 == "transcript", "", gsub('.*exon_number \\"(\\d+)\\".*', "\\1", V9)))]
  updatePhase <- gtf2[V8 %in% c("1", "2") & V3 == "CDS" & exonnumber == 1]
  updatePhase2 <- unique(updatePhase[, .(geneid = .SD$geneid, transcriptid = .SD$transcriptid, exonnumber = c(.SD$exonnumber, .SD$exonnumber, ""), V3 = c("CDS", "exon", "transcript"), ustart = .SD$V4 - 3 + as.numeric(.SD$V8), uend = .SD$V5 + 3 - as.numeric(.SD$V8)), by = seq(1:nrow(updatePhase))])[, .(ustart = min(.SD$ustart), uend = max(.SD$uend)), by = c("geneid", "transcriptid", "exonnumber", "V3")]
  gtf3 <- merge(gtf2, updatePhase2, by = c("geneid", "transcriptid", "exonnumber", "V3"), all.x = T)
  gtf3[, V4p := ifelse(V7 == "+" & !is.na(ustart), ifelse(V4 > ustart, ustart, V4), V4)][, V5p := ifelse(V7 == "-" & !is.na(uend), ifelse(V5 < uend, uend, V5), V5)]
  gtf4 <- gtf3[, .(V1, V2, V3, V4 = as.integer(V4p), V5 = as.integer(V5p), V6, V7, V8, V9)]
  header <- fread("./data/Homo_sapiens.GRCh37.75.gtf", header = F, nrows = 5, sep = ",")
  fwrite(header, file = "./data/Homo_sapiens.GRCh37.75_2.gtf", col.names = F, append = F, eol = "\n", quote = F)
  fwrite(gtf4, file = "./data/Homo_sapiens.GRCh37.75_2.gtf", col.names = F, append = T, sep = "\t", quote = F, eol = "\n")
  Txdb <- makeTxDbFromGFF(file = "./data/Homo_sapiens.GRCh37.75_2.gtf", format = "gtf", dataSource = "Ensembl", organism="Homo sapiens")
  saveDb(Txdb, file = "./data/Txdb.sqlite")
} else Txdb <- AnnotationDbi::loadDb(file = "./data/Txdb.sqlite")
# Get cannonical transcripts: (download from http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh37/appris_data.principal.txt)
cano <- fread("data/hgTables_canonical_transcript.txt", header = F, sep = "\t")
cano2 <- data.table(TXID = as.character(transcripts(Txdb)$tx_id), tx_name = transcripts(Txdb)$tx_name)[, is_cano := ifelse(tx_name %in% gsub("\\..*", "", cano[grepl("PRINCIPAL", V5)]$V3), "Y", "")]

# Function to get access to exac rest api
restAPI <- function(query){
  base <- "http://exac.hms.harvard.edu/"
  call <- paste0(base, query)
  get_res <- GET(call)
  get_res_text <- httr::content(get_res, "text", encoding = "UTF-8")
  res_json <- fromJSON(get_res_text, flatten = T)
  return(res_json)
}

# Function to parse CIGAR string and extract ExAC database friendly variants
get_relevant_var <- function(chr, start, end, ref, alt, width, type, cigar) {
  get_relevant_var_snpindel <- function(chr, start, end, ref, alt, width, type, cigar) {
    if (grepl("^[0-9]+X", cigar)) offset = 0 else offset <- as.numeric(gsub("^([0-9]+)M.*", "\\1", cigar))
    if (!grepl("[ID]", cigar)) insdeln = 0 else insdeln <- as.numeric(gsub("^[0-9]+M([0-9]+)[ID].*", "\\1", cigar))
    if (type == "snp") offset = offset + 1
    starttmp <- start + offset - 1
    reftmp <- substring(ref, offset, offset + ifelse(type == "del", insdeln, 0))
    alttmp <- substring(alt, offset, offset + ifelse(type == "ins", insdeln, 0))
    var_id <- paste(chr, "-", starttmp, "-", reftmp, "-", alttmp, sep = "")
    return(var_id)
  }
  get_relevant_var_mnp <- function(chr, start, end, ref, alt, cigar_expand) {
    starttmp <- seq(start, end, 1)[cigar_expand[cigar_expand != "I"] == "X"]
    reftmp <- strsplit(ref, "")[[1]][cigar_expand[cigar_expand != "I"] == "X"]
    alttmp <- strsplit(alt, "")[[1]][cigar_expand[cigar_expand != "D"] == "X"]
    var_id <- paste(chr, "-", starttmp, "-", reftmp, "-", alttmp, sep = "")
    return(var_id)
  }
  
  if (type %in% c("snp", "ins", "del")) {
    return(get_relevant_var_snpindel(chr, start, end, ref, alt, width, type, cigar))
  }
  if (type == "mnp") {
    cigar_expand <- rep("X", as.numeric(gsub("^([0-9]+)X.*", "\\1", cigar)))
    return(get_relevant_var_mnp(chr, start, end, ref, alt, cigar_expand))
  }
  if (type == "complex") {
    var_id <- NULL
    mat <- data.table(V1 = as.numeric(strsplit(cigar, "[A-Z]+")[[1]]), 
                      V2 = strsplit(cigar, "[0-9]+")[[1]][strsplit(cigar, "[0-9]+")[[1]] != ""]
    )[, V3 := apply(.SD, 1, function(x){rep(x[2], x[1])})]
    if (grepl("[ID]", cigar)) {
      newnm <- nchar(strsplit(paste(gsub("X", "M", unlist(mat$V3)), collapse = ""), "[ID]+")[[1]])[1]
      insdeln <- as.numeric(gsub(".*[0-9]+M([0-9]+)[ID].*", "\\1", cigar))
      cigar <- paste0(newnm, "M", insdeln, ifelse(grepl("I", cigar), "I", "D"))
      var_id <- rbind(var_id, get_relevant_var_snpindel(chr, start, end, ref, alt, width, type, cigar))
    }
    var_id <- rbind(var_id, get_relevant_var_mnp(chr, start, end, ref, alt, unlist(mat$V3)))
    return(var_id)
  }
}


# Function to parse vcf file and perform amino acid annotation
parse_vcf <- function(vcf_file) {
  # Read in vcf file
  withProgress(message = 'Load vcf file...', value = 0, {
    vcf <- readVcf(vcf_file, "hg19")
    incProgress(1/5, detail = paste(1, "/", 5, sep=""))
    # Get position matrix
    pos <- rowRanges(vcf)
    pos_df <- as.data.table(pos)[, ALT := lapply(ALT, function(x){as.character(x)})]
    incProgress(2/5, detail = paste(2, "/", 5, sep=""))
    # Get format matrix
    format <- geno(vcf)
    format_df <- do.call(cbind, lapply(1:length(format), function(x) {
      res <- as.data.table(format[[x]]); colnames(res) <- paste(names(format)[x], "_", colnames(res), sep = "");   
      for(i in 1:ncol(res)) attr(res[[i]], "label") <- geno(header(vcf))$Description[x]; return(res)
    }))
    incProgress(3/5, detail = paste(3, "/", 5, sep=""))
    # Get info matrix
    info <- as.data.table(info(vcf))
    for(x in 1:ncol(info)) {
      attr(info[[x]], "label") <- info(header(vcf))$Description[x]
    }
    incProgress(4/5, detail = paste(4, "/", 5, sep=""))
    
    all_df <- cbind(cbind(pos_df, format_df), info)
    all_df[, no_of_ref := unlist(lapply(ALT, length))]
    incProgress(5/5, detail = paste(5, "/", 5, sep=""))
  }) 
  
  # loc <- as.data.table(locateVariants(vcf, Txdb, AllVariants()))
  # Get coding annotations
  withProgress(message = 'Variants being annotated...', value = 0, {
    coding <- predictCoding(vcf, keepStandardChromosomes(Txdb, "Homo sapiens"), seqSource=BSgenome.Hsapiens.1000genomes.hs37d5, ignore.strand = T)
    incProgress(1/5, detail = paste(1, "/", 5, sep=""))
    coding_df <- merge(as.data.table(coding)[, HGVSp := ifelse(.SD$CONSEQUENCE != "frameshift", paste(paste("p.", paste(c(AMINO_ACID_CODE, "*" = "Ter")[strsplit(.SD$REFAA, "")[[1]]], collapse = ""), paste(.SD$PROTEINLOC[[1]], collapse = "-"), paste(c(AMINO_ACID_CODE, "*" = "Ter")[strsplit(.SD$VARAA, "")[[1]]], collapse = ""), sep = ""), collapse = ";"), ""), by = seq(1:length(coding))], cano2, by = "TXID", all.x = T)[, CONSEQUENCE2 := ifelse(CONSEQUENCE == "nonsynonymous", ifelse(nchar(REFAA) == nchar(VARAA), ifelse(!grepl("Ter", HGVSp), "Substitution", ifelse(grepl("^p.Ter", HGVSp), "Stop Lost", "Stop Gain")), ifelse(nchar(REFAA) < nchar(VARAA), "Insertion", "Deletion")), as.character(CONSEQUENCE))]
    incProgress(2/5, detail = paste(2, "/", 5, sep=""))
    # Derive worst consequence based on the logic: frameshift > nonsense > stop gain > stop lost > deletion > insertion > substitution > synonymous
    coding_df2 <- coding_df[, WCONSEQUENCE := ifelse(sum(grepl("frameshift", .SD$CONSEQUENCE2)) > 0, "frameshift", ifelse(sum(grepl("nonsense", .SD$CONSEQUENCE2)) > 0, "nonsense", ifelse(sum(grepl("Stop Gain", .SD$CONSEQUENCE2)) > 0, "Stop Gain", ifelse(sum(grepl("Stop Lost", .SD$CONSEQUENCE2)) > 0, "Stop Lost", ifelse(sum(grepl("Deletion", .SD$CONSEQUENCE2)) > 0, "Deletion", ifelse(sum(grepl("Insertion", .SD$CONSEQUENCE2)) > 0, "Insertion", ifelse(sum(grepl("Substitution", .SD$CONSEQUENCE2)) > 0, "Substitution", ifelse(sum(grepl("synonymous", .SD$CONSEQUENCE2)) > 0, "synonymous", NA)))))))), by = c("seqnames", "start", "end", "REF", "varAllele")][CONSEQUENCE2 == WCONSEQUENCE][, hasCanonical := ifelse(sum(.SD$is_cano == "Y") > 0, "Y", ""), by = c("seqnames", "start", "end", "REF", "varAllele")]
    incProgress(3/5, detail = paste(3, "/", 5, sep=""))
    coding_df2[, WHGVSp := ifelse(unique(.SD$WCONSEQUENCE) == "frameshift", "", paste(sort(unique(.SD[(hasCanonical == "Y" & is_cano == "Y") | hasCanonical == ""]$HGVSp)), collapse = ";")), by = c("seqnames", "start", "end", "REF", "varAllele")]
    coding_df3 <- coding_df2[, .(varAllele = paste(unique(.SD[, .(varAllele, WCONSEQUENCE, WHGVSp)])$varAllele, collapse = ","), CONSEQUENCE = paste(sort(unique(.SD$CONSEQUENCE2)), collapse = ";"), WCONSEQUENCE = paste(unique(.SD[, .(varAllele, WCONSEQUENCE, WHGVSp)])$WCONSEQUENCE, collapse = ","), WHGVSp = paste(unique(.SD[, .(varAllele, WCONSEQUENCE, WHGVSp)])$WHGVSp, collapse = ",")), by = c("seqnames", "start", "end", "REF")]
    incProgress(4/5, detail = paste(4, "/", 5, sep=""))
    all_df2 <- merge(all_df, coding_df3, by = c("seqnames", "start", "end", "REF"), all = T)
    all_df2[, order := 1:nrow(all_df2)]
    incProgress(5/5, detail = paste(5, "/", 5, sep=""))
  })
  
  # query the Exac database for all the variant information (takes ~7 min to run)
  withProgress(message = 'Query relevant variant information from ExAC database...', value = 0, {
    listres <- NULL
    for (i in 1:nrow(all_df2)) {
      incProgress(1/nrow(all_df2), detail = paste(i, "/", nrow(all_df2), sep=""))
      chr <- all_df2$seqnames[i]
      start <- all_df2$start[i]
      end <- all_df2$end[i]
      res <- restAPI(paste0("/rest/region/variants_in_region/", chr, "-", start, "-", end))
      listres <- c(listres, list(res))
    }
  })
  # filter out for the relevant variant to the best extend
  withProgress(message = 'Filter out unrelevant variants from queried information...', value = 0, {
    listres2 <- NULL
    relev_var <- NULL
    for (i in 1:length(listres)) {
      incProgress(1/length(listres), detail = paste(i, "/", length(listres), sep=""))
      chr <- all_df2$seqnames[i]
      start <- all_df2$start[i]
      end <- all_df2$end[i]
      ref <- all_df2$REF[i]
      alt <- all_df2$ALT[[i]]
      width <- all_df2$width[i]
      type <- all_df2$TYPE[[i]]
      cigar <- all_df2$CIGAR[[i]]
      if (length(listres[[i]]) == 0) {
        listres2 <- c(listres2, list(as.data.table(listres[[i]])))
        relev_var <- c(relev_var, NA)
      } else {
        res <- as.data.table(listres[[i]])
        res_final <- NULL
        var_final <- NULL
        for (j in 1:length(alt)) {
          var_id <- get_relevant_var(chr, start, end, ref, alt[j], width, type[j], cigar[j])
          var_final <- c(var_final, as.character(ifelse(sum(var_id %in% res$variant_id) > 0, paste(var_id[var_id %in% res$variant_id], collapse = ";"), NA)))
          restmp <- res[variant_id %in% var_id & (is.null(res_final) | !variant_id %in% res_final$variant_id)]
          res_final <- rbind(res_final, restmp)
        }
        relev_var <- c(relev_var, list(var_final))
        listres2 <- c(listres2, list(res_final))
      }
    }
  })
  # Combine Exac extracted data and reformat the data table
  withProgress(message = 'Combine Exac extracted data and reformat the data table...', value = 0, {
    all_df2[, Exac := listres2][, ExAC_var := relev_var]
    incProgress(1/6, detail = paste(1, "/", 6, sep=""))
    all_df2[, ':='(allele_freq = ifelse(nrow(.SD$Exac[[1]]) > 0 & "allele_freq" %in% colnames(.SD$Exac[[1]]), paste(sapply(.SD$ExAC_var[[1]], function(x){paste(round(as.data.table(.SD$Exac[[1]])[variant_id %in% strsplit(x, ";")[[1]]]$allele_freq, digits = 3), collapse = ";")}), collapse = ","), "")), by = seq(1:nrow(all_df2))]
    incProgress(2/6, detail = paste(2, "/", 6, sep=""))
    all_df2[, ':='(category = ifelse(nrow(.SD$Exac[[1]]) > 0 & "category" %in% colnames(.SD$Exac[[1]]), paste(sapply(.SD$ExAC_var[[1]], function(x){paste(as.data.table(.SD$Exac[[1]])[variant_id %in% strsplit(x, ";")[[1]]]$category, collapse = ";")}), collapse = ","), "")), by = seq(1:nrow(all_df2))]
    incProgress(3/6, detail = paste(3, "/", 6, sep=""))
    all_df2[, ':='(major_consequence = ifelse(nrow(.SD$Exac[[1]]) > 0 & "major_consequence" %in% colnames(.SD$Exac[[1]]), paste(sapply(.SD$ExAC_var[[1]], function(x){paste(as.data.table(.SD$Exac[[1]])[variant_id %in% strsplit(x, ";")[[1]]]$major_consequence, collapse = ";")}), collapse = ","), "")), by = seq(1:nrow(all_df2))]
    incProgress(4/6, detail = paste(4, "/", 6, sep=""))
    all_df2[, ':='(HGVSp = ifelse(nrow(.SD$Exac[[1]]) > 0 & "HGVSp" %in% colnames(.SD$Exac[[1]]), paste(sapply(.SD$ExAC_var[[1]], function(x){paste(as.data.table(.SD$Exac[[1]])[variant_id %in% strsplit(x, ";")[[1]]]$HGVSp, collapse = ";")}), collapse = ","), "")), by = seq(1:nrow(all_df2))]
    incProgress(5/6, detail = paste(5, "/", 6, sep=""))
    
    # Format for the final output
    final_dat <- all_df2[, .(chr = seqnames, 
                             start, 
                             end, 
                             REF, 
                             ALT = unlist(lapply(ALT, function(x){paste(x, collapse = ",")})), 
                             TYPE = unlist(lapply(TYPE, function(x){paste(x, collapse = ",")})), 
                             CIGAR = unlist(lapply(CIGAR, function(x){paste(x, collapse = ",")})), 
                             WCONSEQUENCE = ifelse(is.na(WCONSEQUENCE), unlist(lapply(TYPE, function(x){paste(rep("Intergenic", length(x)), collapse = ",")})), ifelse(grepl("synonymous", WCONSEQUENCE), gsub("synonymous", "Silent", WCONSEQUENCE), WCONSEQUENCE)), 
                             DP, 
                             AO = unlist(lapply(AO, function(x){paste(x, collapse = ",")})), 
                             AF = unlist(lapply(AF, function(x){paste(x, collapse = ",")})), 
                             ExAC_var = unlist(lapply(ExAC_var, function(x){paste(x, collapse = ",")})), 
                             allele_freq, 
                             category, 
                             major_consequence, 
                             HGVSp)]
    incProgress(6/6, detail = paste(6, "/", 6, sep=""))
  }) 
  return(final_dat)
}