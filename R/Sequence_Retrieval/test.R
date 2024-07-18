

species_name <- "Oryza sativa"


a1 <- entrez_search(db = "taxonomy",
              term = species_name) %>% 
  {entrez_link(dbfrom = "taxonomy",
              db = "assembly",
              id =. $ids)} %>% 
  {.$links$taxonomy_assembly} %>%
  split(., ceiling(seq_along(.) / 100)) %>% 
  map(~ entrez_summary(db = "assembly", id = .)) %>% 
  flatten() %>% 
  map_df(~ tibble(
    Assembly.UID = coalesce(extract_from_esummary(., "uid"), NA),
    Assembly.Accession = coalesce(extract_from_esummary(., "assemblyaccession"), NA),
    Assembly.Name = coalesce(extract_from_esummary(., "assemblyname"), NA),
    Assembly.Type = coalesce(extract_from_esummary(., "assemblytype"),  NA),
    Assembly.Status = coalesce(extract_from_esummary(., "assemblystatus"), NA),
    Reference.Seq = coalesce(extract_from_esummary(., "refseq_category"), NA),
    Last.Update.Date = coalesce(extract_from_esummary(., "lastupdatedate"), NA),
    Submitter = coalesce(extract_from_esummary(., "submitterorganization"), NA),
    Bioprojects = coalesce(extract_from_esummary(., "gb_bioprojects")[1], NA),
    Biosample = coalesce(extract_from_esummary(., "biosampleaccn"), NA),
    Organism = coalesce(extract_from_esummary(., "organism"), NA),
    URL = coalesce(extract_from_esummary(., "ftppath_genbank"), NA)
  ))


assembly_accession <- c("GCA_001580425.1")

b <- entrez_search(db = "assembly",
                   term = assembly_accession) %>% 
  {entrez_link(dbfrom = "assembly",
               db = "nuccore",
               id =. $ids)} %>% 
  {.$links$assembly_nuccore} %>%
  split(., ceiling(seq_along(.) / 100)) %>% 
  map(~ entrez_summary(db = "nuccore", id = .)) %>% 
  flatten() %>% 
  map_df(~ tibble(
    Sequence.UID = coalesce(extract_from_esummary(., "uid"), NA),
    Sequence.Accession = coalesce(extract_from_esummary(., "accessionversion"), NA),
    Sequence.Source = coalesce(extract_from_esummary(., "sourcedb"), NA),
    Sequence.Name = coalesce(extract_from_esummary(., "title"), NA),
    Sequence.Type = coalesce(extract_from_esummary(., "biomol"),  NA),
    Sequence.Structure = coalesce(extract_from_esummary(., "genome"), NA),
    Sequence.Completeness = coalesce(extract_from_esummary(., "completeness"), NA),
    Sequence.Length = coalesce(extract_from_esummary(., "slen"), NA),
    Last.Update.Date = coalesce(extract_from_esummary(., "updatedate"), NA),
    Biosample = coalesce(extract_from_esummary(., "biosample"), NA),
    Organism = coalesce(extract_from_esummary(., "organism"), NA),
    Strain = coalesce(extract_from_esummary(., "strain"), NA)) 
  )


c <- entrez_search(db = "taxonomy",
                   term = species_name) %>% 
  {entrez_link(dbfrom = "taxonomy",
               db = "genome",
               id =. $ids)} %>% 
  {entrez_link(dbfrom = "genome",
               db = "nuccore",
               id =. $links$taxonomy_genome)} %>% 
  {.$links$genome_nuccore} %>%
  split(., ceiling(seq_along(.) / 100)) %>% 
  map(~ entrez_summary(db = "nuccore", id = ., retmax = 10)) %>% 
  flatten() %>% 
  map_df(~ tibble(
    Sequence.UID = coalesce(extract_from_esummary(., "uid"), NA),
    Sequence.Accession = coalesce(extract_from_esummary(., "accessionversion"), NA),
    Sequence.Source = coalesce(extract_from_esummary(., "sourcedb"), NA)
  ))


