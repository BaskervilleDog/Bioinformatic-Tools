

species_name <- "Mus musculus"
file_path <- "refseq_files"
refseq_filter <- 'AND "latest refseq"[filter])'
genbank_filter <- 'AND ("latest genbank"[filter] NOT "latest refseq"[filter])'




assembly_explorer <- function(species_name, database_type = c("refseq", "genbank"), retmax = 100) {
  database_type <- match.arg(database_type)
  search_term <- paste0('"', species_name, '"', '[Organism]')
  
  if (database_type == "refseq") {
    search_term <- paste0(search_term, ' AND "latest refseq"[filter]')
  } else if (database_type == "genbank") {
    search_term <- paste0(search_term, 'AND ("latest genbank"[filter] NOT "latest refseq"[filter])')
  }
  
  search_results <- entrez_search(db = "assembly", term = search_term, retmax = retmax)
  
  if (length(search_results$ids) == 0) {
    message("No results found for the specified species and database type.")
    return(NULL)
  }
  
  search_results$ids %>% 
    split(., ceiling(seq_along(.) / 100)) %>% 
    map(~ if(length(.) == 1) list(entrez_summary(db = "assembly", id = .)) else entrez_summary(db = "assembly", id = .)) %>%   
    flatten() %>% 
    map_df(~ tibble(
      Assembly.UID = coalesce(extract_from_esummary(., "uid"), NA),
      Assembly.Accession = coalesce(extract_from_esummary(., "assemblyaccession"), NA),
      Assembly.Name = coalesce(extract_from_esummary(., "assemblyname"), NA),
      Assembly.Type = coalesce(extract_from_esummary(., "assemblytype"),  NA),
      Assembly.Status = coalesce(extract_from_esummary(., "assemblystatus"), NA),
      Reference.Status = coalesce(extract_from_esummary(., "refseq_category"), NA),
      Last.Update.Date = coalesce(extract_from_esummary(., "lastupdatedate"), NA),
      Submitter = coalesce(extract_from_esummary(., "submitterorganization"), NA),
      Bioprojects = coalesce(extract_from_esummary(., "gb_bioprojects")$bioprojectaccn, NA),
      Biosample = coalesce(extract_from_esummary(., "biosampleaccn"), NA),
      Organism = coalesce(extract_from_esummary(., "organism"), NA),
      Taxonomy.Id = coalesce(extract_from_esummary(., "taxid"), NA),
      Species.Name = coalesce(extract_from_esummary(., "speciesname"), NA),
      URL.GenBank = coalesce(extract_from_esummary(., "ftppath_genbank"), NA),
      URL.RefSeq = coalesce(extract_from_esummary(., "ftppath_refseq"), NA)
    ))
}

assembly_explorer("Oryza sativa", "refseq", retmax = 200)



get_refseq <- function(species_name,
                       database_type = c("refseq", "genbank"),
                       type_of_assembly = c("Complete Genome","Chromosome", "Scaffold", "Contig"),
                       retmax = 100,
                       file_path) {
  database_type <- match.arg(database_type)
  search_term <- paste0('"', species_name, '"', '[Organism]')
  
  if (database_type == "refseq") {
    search_term <- paste0(search_term, ' AND "latest refseq"[filter]')
  } else if (database_type == "genbank") {
    search_term <- paste0(search_term, ' AND ("latest genbank"[filter] NOT "latest refseq"[filter])')
  }
  
  if (!dir.exists(file_path)) {
    dir.create(file_path, recursive = TRUE)
  }
  
  tryCatch({
    URL <- entrez_search(db = "assembly", term = search_term, retmax = retmax)
    
    if (length(URL$ids) == 0) {
      message("No results found for the specified species and database type.")
      return(NULL)
    }
    
    if (database_type == "refseq") {
      URL$ids %>% 
        split(., ceiling(seq_along(.) / 100)) %>% 
        map(~ if(length(.) == 1) list(entrez_summary(db = "assembly", id = .)) else entrez_summary(db = "assembly", id = .)) %>%   
        flatten() %>% 
        map_df(~ tibble(
          Assembly.Accession = coalesce(extract_from_esummary(., "assemblyaccession"), NA),
          URL.RefSeq = paste0(coalesce(extract_from_esummary(., "ftppath_refseq"), NA), "/")
        )) -> URL
      
    } else if (database_type == "genbank"){
      URL$ids %>% 
        split(., ceiling(seq_along(.) / 100)) %>% 
        map(~ if(length(.) == 1) list(entrez_summary(db = "assembly", id = .)) else entrez_summary(db = "assembly", id = .)) %>%   
        flatten() %>% 
        map_df(~ tibble(
          Assembly.Accession = coalesce(extract_from_esummary(., "assemblyaccession"), NA),
          Assembly.Status = coalesce(extract_from_esummary(., "assemblystatus"), NA),
          URL.RefSeq = paste0(coalesce(extract_from_esummary(., "ftppath_genbank"), NA), "/")
        )) %>% 
        filter(Assembly.Status %in% type_of_assembly) -> URL
    }
    
    if (nrow(URL) == 0) {
      stop("No Assemblies were found for the given species.")
    }
    
    # Generate and download file lists for each URL.RefSeq
    lapply(seq_len(nrow(URL)), function(i) {
      refseq_url <- URL$URL.RefSeq[i]
      accession <- URL$Assembly.Accession[i]
      subdir_path <- file.path(file_path, accession)
      
      if (!dir.exists(subdir_path)) {
        dir.create(subdir_path, recursive = TRUE)
      }
      
      file_list <- getURL(refseq_url, ftp.use.epsv = F, dirlistonly = T) %>%
        {strsplit(., "\n")[[1]]} %>%
        gsub("\r", "", .) %>%
        list(
          meta = grep("report.txt|regions.txt|stats.txt|status.txt", ., value = TRUE),
          genomic = grep("genomic.fna|genomic.gff|genomic.gtf|genomic.gbff|_gaps|cds_from_genomic", ., value = TRUE),
          rna = grep("rna.fna|rna_from_genomic.fna|rna_gbff", ., value = TRUE),
          protein = grep("_protein", ., value = TRUE)
        )
      
      for (category in names(file_list)) {
        for (file in file_list[[category]]) {
          download.file(url = paste0(refseq_url, file), destfile = file.path(subdir_path, file))
        }
      }
    })
    
  }, error = function(e) {
    message("An error occurred during processing: ", e$message)
  })
}

get_refseq("Saccharomyces cerevisiae", "refseq", file_path = "Assembly_files")






aa2 <- entrez_search(db = "assembly", term = paste(paste0('"',species_name,'"', '[Organism]'), genbank_filter), retmax = 100) %>% 
  {.$ids} %>% 
  split(., ceiling(seq_along(.) / 100)) %>% 
  map(~ if(length(.) == 1) list(entrez_summary(db = "assembly", id = .)) else entrez_summary(db = "assembly", id = .)) %>%   
  flatten() %>% 
  map_df(~ tibble(
    Assembly.UID = coalesce(extract_from_esummary(., "uid"), NA),
    Assembly.Accession = coalesce(extract_from_esummary(., "assemblyaccession"), NA),
    Assembly.Name = coalesce(extract_from_esummary(., "assemblyname"), NA),
    Assembly.Type = coalesce(extract_from_esummary(., "assemblytype"),  NA),
    Assembly.Status = coalesce(extract_from_esummary(., "assemblystatus"), NA),
    Reference.Status = coalesce(extract_from_esummary(., "refseq_category"), NA),
    Last.Update.Date = coalesce(extract_from_esummary(., "lastupdatedate"), NA),
    Submitter = coalesce(extract_from_esummary(., "submitterorganization"), NA),
    Bioprojects = coalesce(extract_from_esummary(., "gb_bioprojects")[1], NA),
    Biosample = coalesce(extract_from_esummary(., "biosampleaccn"), NA),
    Organism = coalesce(extract_from_esummary(., "organism"), NA),
    Taxonomy.Id = coalesce(extract_from_esummary(., "taxid"), NA),
    Species.Name = coalesce(extract_from_esummary(., "speciesname"), NA),
    URL.GenBank = coalesce(extract_from_esummary(., "ftppath_genbank"), NA),
    URL.RefSeq = coalesce(extract_from_esummary(., "ftppath_refseq"), NA),
    URL.Annotation = coalesce(extract_from_esummary(., "annotrpturl"), NA),
    URL.Assembly.Statistics = coalesce(extract_from_esummary(., "ftppath_stats_rpt"), NA),
    URL.Assembly.Report = coalesce(extract_from_esummary(., "ftppath_assembly_rpt"), NA)
  ))


assembly_name <- "UU_Cfam_GSD_1.0"


genecount <- entrez_search(db = "gene", term = paste(paste0(assembly_name, "[All fields]"), "AND (source_genomic[properties] AND srcdb refseq[Properties]"), retmax = 999999) %>% 
  {.$ids} %>% 
  length(.)


a2 <- entrez_search(db = "gene", term = paste(paste0(assembly_name, "[All fields]"), "AND (source_genomic[properties] AND srcdb refseq[Properties]"), retmax = 1000) %>% 
  {.$ids} %>% 
  split(., ceiling(seq_along(.) / 100)) %>% 
  map(~ if(length(.) == 1) list(entrez_summary(db = "gene", id = .)) else entrez_summary(db = "gene", id = .)) %>%  
  flatten() %>% 
  map_df(~ tibble(
    Gene.UID = coalesce(extract_from_esummary(., "uid"), NA),
    Gene.Accession = coalesce(extract_from_esummary(., "genomicinfo")$chraccver, NA),
    Gene.Name = coalesce(extract_from_esummary(., "name"), NA),
    Gene.Chromosome = coalesce(extract_from_esummary(., "chromosome"), NA),
    Gene.Start = coalesce(extract_from_esummary(., "genomicinfo")$chrstart, NA),
    Gene.Stop = coalesce(extract_from_esummary(., "genomicinfo")$chrstop, NA),
    Gene.Description = coalesce(extract_from_esummary(., "description"), NA),
    Gene.Organism = coalesce(extract_from_esummary(., "organism")$scientificname, NA),
    Gene.Taxonomy = coalesce(extract_from_esummary(., "organism")$taxid, NA)
  ))


a3 <- entrez_search(db = "gene", term = paste(paste0(assembly_name, "[All fields]"), "AND (source_genomic[properties] AND srcdb refseq[Properties]"), retmax = 1000) %>% 
  {.$ids} %>% 
  split(., ceiling(seq_along(.) / 100)) %>% 
  map(~ if(length(.) == 1) list(entrez_summary(db = "gene", id = .)) else entrez_summary(db = "gene", id = .)) %>%  
  flatten() %>% 
  map_df(~ tibble(
    Gene.UID = coalesce(extract_from_esummary(., "uid"), NA),
    Gene.Accession = coalesce(extract_from_esummary(., "genomicinfo")$chraccver, NA),
    Gene.Name = coalesce(extract_from_esummary(., "name"), NA),
    Gene.Chromosome = coalesce(extract_from_esummary(., "chromosome"), NA),
    Gene.Start = coalesce(extract_from_esummary(., "genomicinfo")$chrstart, NA),
    Gene.Stop = coalesce(extract_from_esummary(., "genomicinfo")$chrstop, NA),
    Gene.Description = coalesce(extract_from_esummary(., "description"), NA),
    Gene.Organism = coalesce(extract_from_esummary(., "organism")$scientificname, NA),
    Gene.Taxonomy = coalesce(extract_from_esummary(., "organism")$taxid, NA)
  ))
