# Funtions for biological sequence retrieval
# Author: Gianlucca de Urzêda Alves


#==============================================================================#
#==============================================================================#

# Install and load "rentrez" package (NCBI EUtils R package)
if (!requireNamespace("rentrez", quietly = TRUE)) {
  install.packages("rentrez")
}

library(rentrez)

# Install and load "tidyverse" package (Data manipulation packages)
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

library(tidyverse)

# Install and load "biomaRt" package (ENSEMBL Data R package in BiocondutR)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

library(biomaRt)

if (!requireNamespace("BSgenome", quietly = TRUE))
  BiocManager::install("BSgenome")

library(BSgenome)

if (!requireNamespace("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")

library(Biostrings)

if (!requireNamespace("GenomicRanges", quietly = TRUE))
  BiocManager::install("GenomicRanges")

library(GenomicRanges)

if (!requireNamespace("RCurl", quietly = TRUE))
  install.packages("RCurl")

library(RCurl)

#==============================================================================#
#==============================================================================#

# LiteratureRetrieval allows the user to obtain a tibble containing the main
# aspects of information regarding Scientific Literature

LiteratureRetrieval <- function(term, db = "pubmed", retmax = 10000) {
  # Validate inputs
  if (!is.character(term) || length(term) != 1) {
    stop("The 'term' parameter must be a single string.")
  }
  if (!is.character(db) || length(db) != 1) {
    stop("The 'db' parameter must be a single string.")
  }
  if (!is.numeric(retmax) || retmax <= 0) {
    stop("The 'retmax' parameter must be a positive number.")
  }
  
  # Retrieve literature results
  literature_results <- entrez_search(db = db, term = term, retmax = retmax, use_history = TRUE)
  
  # Check if any literature was found
  if (length(literature_results$ids) == 0) {
    stop("No literature found for the specified search terms.")
  }
  
  # Function to retrieve summaries in batches using web history
  get_summary_batch <- function(web_history, retstart, retmax) {
    literature_summary <- entrez_summary(db = db, web_history = web_history, retstart = retstart, retmax = retmax)
    tibble(
      Uid = extract_from_esummary(literature_summary, "uid", simplify = TRUE),
      Publication.Date = str_extract(extract_from_esummary(literature_summary, "pubdate", simplify = TRUE), "\\d{4}"),
      Journal = extract_from_esummary(literature_summary, "fulljournalname", simplify = TRUE),
      ISSN = extract_from_esummary(literature_summary, "issn", simplify = TRUE),
      ESSN = extract_from_esummary(literature_summary, "essn", simplify = TRUE),
      Title = extract_from_esummary(literature_summary, "title", simplify = TRUE),
      Electronic.Location = extract_from_esummary(literature_summary, "elocationid", simplify = TRUE),
      Authors = sapply(literature_summary, function(x) paste(x$authors$name, collapse = ", ")),
      Publication.Type = sapply(literature_summary, function(x) paste(x$pubtype, collapse = ", "))
    )
  }
  
  # Determine the number of batches needed
  num_batches <- ceiling(length(literature_results$ids) / 500)
  
  # Retrieve summaries for each batch and combine into a single tibble
  literature_table <- do.call(rbind, lapply(0:(num_batches - 1), function(i) {
    get_summary_batch(literature_results$web_history, retstart = i * 500, retmax = 500)
  }))
  
  return(literature_table)
}

# Example of usage:
# Specify a search term such as "Cellulose synthase in Arabidopsis thaliana"
# Specify a NCBI database for Scientific Literature ("pubmed"); default = "pubmed"
# Specify a number of hits to display; default = 10000
# Lit <- LiteratureRetrieval("Cellulose synthase in Arabidopsis thaliana", db = "pubmed", retmax = 250)

# Lit <- LiteratureRetrieval("Cellulose synthase in Arabidopsis thaliana") # default search parameters

#==============================================================================#
#==============================================================================#

# GenomeRetrieval allows the user to obtain a tibble containing the main
# aspects of information regarding Genomes of a specified organism

GenomeRetrieval <- function(organism_name, retmax = 10000) {
  # Search for the taxonomy ID of the given organism
  taxonomy <- entrez_search(db = "taxonomy", term = organism_name, retmax = 1)
  taxonomy_id <- taxonomy$ids
  
  if (length(taxonomy$ids) == 0) {
    stop("No Taxonomy ID found for the specified organism.")
  }
  
  # Link the taxonomy ID to the genome database
  link_taxonomy_genome <- entrez_link(dbfrom = "taxonomy",
                                      id = taxonomy_id,
                                      db = "genome")
  
  genome_ids <- link_taxonomy_genome$links$taxonomy_genome
  
  # Link the genome IDs to the nucleotide database and retrieve web history
  link_genome_sequences <- entrez_link(dbfrom = "genome",
                                       id = genome_ids,
                                       db = "nuccore",
                                       cmd = "neighbor_history")
  
  web_history <- link_genome_sequences$web_histories
  
  # Function to fetch summaries in batches
  fetch_summaries_in_batches <- function(web_history, batch_size, retmax) {
    all_summaries <- list()
    for (start in seq(0, retmax - 1, by = batch_size)) {
      batch_summary <- entrez_summary(db = "nuccore",
                                      web_history = web_history$genome_nuccore,
                                      retstart = start,
                                      retmax = batch_size)
      all_summaries <- c(all_summaries, batch_summary)
      if (length(batch_summary) < batch_size) break
    }
    return(all_summaries)
  }
  
  # Fetch summaries in batches
  batch_size <- 500
  genome_summary <- fetch_summaries_in_batches(web_history, batch_size, retmax)
  
  # Extract function to handle lists of summaries
  extract_from_esummary_list <- function(esummary_list, field, simplify = FALSE) {
    sapply(esummary_list, function(record) extract_from_esummary(record, field, simplify))
  }
  
  # Create a table with the relevant details
  genome_table <- tibble(
    GENOME_Uid = extract_from_esummary_list(genome_summary, "uid", simplify = TRUE),
    Date = extract_from_esummary_list(genome_summary, "createdate", simplify = TRUE),
    Accession = extract_from_esummary_list(genome_summary, "caption", simplify = TRUE),
    Acc_version = extract_from_esummary_list(genome_summary, "accessionversion", simplify = TRUE),
    Title = extract_from_esummary_list(genome_summary, "title", simplify = TRUE),
    Seq.Length = extract_from_esummary_list(genome_summary, "slen", simplify = TRUE),
    Biomolecule = extract_from_esummary_list(genome_summary, "biomol", simplify = TRUE),
    Database = extract_from_esummary_list(genome_summary, "sourcedb", simplify = TRUE),
    Organism = extract_from_esummary_list(genome_summary, "organism", simplify = TRUE)
  )
  
  return(genome_table)
}

# Example of usage:
# Specify a organism common or scientific name such as "Arabidopsis thaliana" or "Rice"
# Specify a number of hits to display; default = 10000
# genome <- GenomeRetrieval("Arabidopsis thaliana", retmax = 10000)

# genome <- GenomeRetrieval("Arabidopsis thaliana") # default search parameters



#==============================================================================#
#==============================================================================#

# GenomeFetcher allows the user to obtain a FASTA file containing the sequences
# of either the Genome UIDs gathered with GenomeRetrieval or by using Genome UIDs from
# NCBI
# **IMPORTANT**: The batch size must be smaller than 50, or the NCBI platform might
# lock you out from using its databases for a few days

GenomeFetcher <- function(genome_ids, fasta_file, batch_size = 40) {
  # Initialize an empty vector to store the fetched sequences
  all_sequences <- character(0)
  
  if (batch_size > 50) {
    stop("Batch size must be 50 or smaller")
  }
  
  # Loop through the genome IDs in batches
  for (i in seq(1, length(genome_ids), by = batch_size)) {
    # Determine the range of the current batch
    batch_end <- min(i + batch_size - 1, length(genome_ids))
    current_batch <- genome_ids[i:batch_end]
    
    # Fetch genome sequences for the current batch
    genome_sequence_fetch <- entrez_fetch(db = "nuccore", id = current_batch, rettype = "fasta")
    
    # Append the fetched sequences to the all_sequences vector
    all_sequences <- c(all_sequences, genome_sequence_fetch)
  }
  
  # Write all fetched sequences to the specified fasta file
  writeLines(all_sequences, fasta_file)
  
  return(paste("Sequences written into", fasta_file))
}

# Example of usage:
# Using genome UIDs gathered with GenomeRetrieval
# The sequences can be written in any format in the final file (FASTA format is recommended)
# batch_size must be smaller than 50
# GenomeFetcher(genome$Uid, "Seq.fasta", batch_size = 40)


#==============================================================================#
#==============================================================================#

# GeneRetrieval allows the user to obtain a tibble containing the main
# Aspects of information regarding specific Genes of a certain organism

GeneRetrieval <- function(gene_names, organism_name, retmax = 2000, batch_size = 500) {
  # Helper function to retrieve gene information for a single gene name
  retrieve_gene_info <- function(gene_name) {
    # Construct the search term
    search_term <- paste0(gene_name, "[TITL] gene AND ", organism_name, "[ORGN]")
    
    # Search for genes in the specified organism
    genes <- entrez_search(db = "gene", term = search_term, retmax = retmax, use_history = TRUE)
    
    # Check if any genes were found
    if (genes$count == 0) {
      warning(paste("No genes found for the gene name:", gene_name))
      return(NULL)
    }
    
    # Split the web history into batches
    web_hist <- genes$web_history
    batches <- seq(0, genes$count - 1, by = batch_size)
    
    # Function to fetch and process summaries for a batch
    fetch_summaries <- function(start) {
      summaries <- entrez_summary(db = "gene", web_history = web_hist, retstart = start, retmax = batch_size)
      
      tibble(
        GeneName = gene_name,
        Uid = extract_from_esummary(summaries, "uid", simplify = TRUE),
        Chr.Acc.Number = sapply(summaries, function(x) x$genomicinfo$chraccver),
        Name = extract_from_esummary(summaries, "name", simplify = TRUE),
        Description = extract_from_esummary(summaries, "description", simplify = TRUE),
        Organism = sapply(summaries, function(x) x$organism$scientificname),
        Chromosome = sapply(summaries, function(x) x$genomicinfo$chrloc),
        Start = sapply(summaries, function(x) x$genomicinfo$chrstart),
        Stop = sapply(summaries, function(x) x$genomicinfo$chrstop),
        Summary = extract_from_esummary(summaries, "summary", simplify = TRUE)
      )
    }
    
    # Fetch and process all batches
    all_genes <- map_dfr(batches, fetch_summaries)
    
    all_genes <- all_genes %>%
      unnest(cols = c(Chr.Acc.Number, Chromosome, Start, Stop))
    
    return(all_genes)
  }
  
  # Iterate over gene names and retrieve information for each
  all_results <- map_dfr(gene_names, retrieve_gene_info)
  
  return(all_results)
}

# Example of usage:
# Specify a gene or list of genes such as "MAPK" or c("MAPK", "BRCA1") and a organism "Homo sapiens"
# Specify a number of hits to display; default = 2000
# The max batch_size for consultation is 500


#genes <- GeneRetrieval("MAPK", "Homo sapiens", retmax = 2000, batch_size = 500)

# genes <- GeneRetrieval(c("MAPK",
#"Alcohol Dehydrogenase",
#"p53",
#"BRCA1",
#"TNF",
#"EGFR",
#"IL6",
#"TGFB1",
#"ESR1",
#"APOE",
#"Transferase"), "Homo sapiens")


#==============================================================================#
#==============================================================================#


#==============================================================================#
#==============================================================================#

#mRNARetrieval <- function(gene_uids) {
  # Link gene IDs to the nuccore database
#  search_links <- entrez_link(dbfrom = "gene", id = gene_uids, db = "nuccore", cmd = "neighbor_history")
  
  # Check if any RNA sequences were found
#  if (length(search_links$web_histories$gene_nuccore_refseqrna) == 0) {
#    stop("No RNA sequences found for the specified gene UIDs.")
#  }
  
  # Fetch summaries of the RNA sequences
#  RNA_summary <- entrez_summary(db = "nuccore", web_history = search_links$web_histories$gene_nuccore_refseqrna)
  
  # Extract details from the summaries
#  RNA_table <- tibble(
#    Uid = extract_from_esummary(RNA_summary, "uid", simplify = TRUE),
#    Accession = extract_from_esummary(RNA_summary, "caption", simplify = TRUE),
#    Acc_version = extract_from_esummary(RNA_summary, "accessionversion", simplify = TRUE),
#    Title = extract_from_esummary(RNA_summary, "title", simplify = TRUE),
#    Seq_Length = extract_from_esummary(RNA_summary, "slen", simplify = TRUE),
#    Biomolecule = extract_from_esummary(RNA_summary, "biomol", simplify = TRUE),
#    Database = extract_from_esummary(RNA_summary, "sourcedb", simplify = TRUE),
#    Organism = extract_from_esummary(RNA_summary, "organism", simplify = TRUE)
#  )
  
  # Return the RNA details table
#  return(RNA_table)
#}


# mRNARetrieval allows the user to obtain a tibble containing the main
# Aspects of information regarding mRNA of previously selected Genes of 
# a certain organism

mRNARetrieval <- function(gene_uids) {
  all_RNA_tables <- list()
  
  for (gene_uid in gene_uids) {
    # Link gene ID to the nuccore database
    search_links <- tryCatch(
      entrez_link(dbfrom = "gene", id = gene_uid, db = "nuccore", cmd = "neighbor_history"),
      error = function(e) {
        warning(paste("Error linking gene ID", gene_uid, ":", e$message))
        return(NULL)
      }
    )
    
    # Check if linking was successful
    if (is.null(search_links)) next
    
    # Check if any RNA sequences were found
    if (length(search_links$web_histories$gene_nuccore_refseqrna) == 0) {
      warning(paste("No RNA sequences found for gene UID", gene_uid))
      next
    }
    
    # Fetch summaries of the RNA sequences
    RNA_summary <- tryCatch(
      entrez_summary(db = "nuccore", web_history = search_links$web_histories$gene_nuccore_refseqrna),
      error = function(e) {
        warning(paste("Error fetching RNA summaries for gene UID", gene_uid, ":", e$message))
        return(NULL)
      }
    )
    
    # Check if fetching summaries was successful
    if (is.null(RNA_summary)) next
    
    # Extract details from the summaries
    RNA_table <- tibble(
      Gene_UID = gene_uid,  # Add the gene UID
      RNA_UID = extract_from_esummary(RNA_summary, "uid", simplify = TRUE),
      Accession = extract_from_esummary(RNA_summary, "caption", simplify = TRUE),
      Acc_version = extract_from_esummary(RNA_summary, "accessionversion", simplify = TRUE),
      Title = extract_from_esummary(RNA_summary, "title", simplify = TRUE),
      Seq_Length = extract_from_esummary(RNA_summary, "slen", simplify = TRUE),
      Biomolecule = extract_from_esummary(RNA_summary, "biomol", simplify = TRUE),
      Database = extract_from_esummary(RNA_summary, "sourcedb", simplify = TRUE),
      Organism = extract_from_esummary(RNA_summary, "organism", simplify = TRUE)
    )
    
    # Append the current RNA table to the list
    all_RNA_tables <- append(all_RNA_tables, list(RNA_table))
  }
  
  # Combine all the RNA tables into one
  combined_RNA_table <- bind_rows(all_RNA_tables)
  
  return(combined_RNA_table)
}

# Example of usage:
# Specify a gene UID or list of gene UIDs (GeneRetrieval might me used for this)
# mRNA <- mRNARetrieval(genes$Uid)

# Or use a specific gene UID
#mRNA <- mRNARetrieval("5594")

# For a single gene, multiple mRNA may exist, therefore the function is expected to retrieve
# an unexpectedly high number of sequences, given a number of gene UIDs. Use a curated and filtered
# gene list for this function

#==============================================================================#
#==============================================================================#

# mRNAFetcher allows the user to obtain a FASTA file containing the sequences
# of either the Gene UIDs gathered with GeneRetrieval or by using Gene UIDs from
# NCBI
# **IMPORTANT**: The batch size must be smaller than 50, or the NCBI platform might
# lock you out from using its databases for a few days

mRNAFetcher <- function(mRNA_ids, fasta_file, batch_size = 40) {
  # Initialize an empty vector to store the fetched sequences
  all_sequences <- character(0)
  
  if (batch_size > 50) {
    stop("Batch size must be 50 or smaller")
  }
  
  # Loop through the mRNA IDs in batches
  for (i in seq(1, length(mRNA_ids), by = batch_size)) {
    # Determine the range of the current batch
    batch_end <- min(i + batch_size - 1, length(mRNA_ids))
    current_batch <- mRNA_ids[i:batch_end]
    
    # Fetch mRNA sequences for the current batch
    mRNA_sequence_fetch <- entrez_fetch(db = "nuccore", id = current_batch, rettype = "fasta")
    
    # Append the fetched sequences to the all_sequences vector
    all_sequences <- c(all_sequences, mRNA_sequence_fetch)
  }
  
  # Write all fetched sequences to the specified fasta file
  writeLines(all_sequences, fasta_file)
  
  return(paste("mRNA Sequences written into", fasta_file))
}

# Example of usage:
# Using gene UIDs gathered with GeneRetrieval
# The sequences can be written in any format in the final file (FASTA format is recommended)
# batch_size must be smaller than 50
# mRNAFetcher(mRNA$RNA_UID,"mRNA_seq.fasta", batch_size = 40)
# mRNAFetcher("1777376004", "mRNA_seq.fasta", batch_size = 40)

#==============================================================================#
#==============================================================================#

# ProteinRetrieval allows the user to obtain a tibble containing the main
# Aspects of information regarding proteins of previously selected Genes of 
# a certain organism

ProteinRetrieval <- function(gene_uids) {
  all_Protein_tables <- list()
  
  for (gene_uid in gene_uids) {
    # Link gene ID to the protein database
    search_links <- tryCatch(
      entrez_link(dbfrom = "gene", id = gene_uid, db = "protein", cmd = "neighbor_history"),
      error = function(e) {
        warning(paste("Error linking gene ID", gene_uid, ":", e$message))
        return(NULL)
      }
    )
    
    # Check if linking was successful
    if (is.null(search_links)) next
    
    # Check if any protein sequences were found
    if (length(search_links$web_histories$gene_protein_refseq) == 0) {
      warning(paste("No protein sequences found for gene UID", gene_uid))
      next
    }
    
    # Fetch summaries of the protein sequences
    Protein_summary <- tryCatch(
      entrez_summary(db = "protein", web_history = search_links$web_histories$gene_protein_refseq),
      error = function(e) {
        warning(paste("Error fetching protein summaries for gene UID", gene_uid, ":", e$message))
        return(NULL)
      }
    )
    
    # Check if fetching summaries was successful
    if (is.null(Protein_summary)) next
    
    # Extract details from the summaries
    Protein_tibble <- tibble(
      Gene_UID = gene_uid,  # Add the gene UID
      PROT_UID = extract_from_esummary(Protein_summary, "uid", simplify = TRUE),
      Accession = extract_from_esummary(Protein_summary, "caption", simplify = TRUE),
      Title = extract_from_esummary(Protein_summary, "title", simplify = TRUE),
      Seq_Length = extract_from_esummary(Protein_summary, "slen", simplify = TRUE),
      Biomolecule = extract_from_esummary(Protein_summary, "moltype", simplify = TRUE),
      Database = extract_from_esummary(Protein_summary, "sourcedb", simplify = TRUE),
      Organism = extract_from_esummary(Protein_summary, "organism", simplify = TRUE)
    )
    
    # Append the current protein table to the list
    all_Protein_tables <- append(all_Protein_tables, list(Protein_tibble))
  }
  
  # Combine all the protein tables into one
  combined_Protein_table <- bind_rows(all_Protein_tables)
  
  return(combined_Protein_table)
}

# Example of usage:
# Specify a gene UID or list of gene UIDs (GeneRetrieval might me used for this)
# protein <- ProteinRetrieval(genes$Uid)

# Or use a specific gene UID
#protein <- ProteinRetrieval("5594")

# For a single gene, multiple proteins may exist, therefore the function is expected to retrieve
# an unexpectedly high number of sequences, given a number of gene UIDs. Use a curated and filtered
# gene list for this function

#====================================================================================================================#
#====================================================================================================================#

# ProteinFetcher allows the user to obtain a FASTA file containing the sequences
# of either the Gene UIDs gathered with GeneRetrieval or by using Gene UIDs from
# NCBI
# **IMPORTANT**: The batch size must be smaller than 50, or the NCBI platform might
# lock you out from using its databases for a few days

ProteinFetcher <- function(protein_ids, fasta_file, batch_size = 40) {
  # Initialize an empty vector to store the fetched sequences
  all_sequences <- character(0)
  
  if (batch_size > 50) {
    stop("Batch size must be 50 or smaller")
  }
  
  # Loop through the protein IDs in batches
  for (i in seq(1, length(protein_ids), by = batch_size)) {
    # Determine the range of the current batch
    batch_end <- min(i + batch_size - 1, length(protein_ids))
    current_batch <- protein_ids[i:batch_end]
    
    # Fetch protein sequences for the current batch
    Protein_sequence_fetch <- entrez_fetch(db = "protein", id = current_batch, rettype = "fasta")
    
    # Append the fetched sequences to the all_sequences vector
    all_sequences <- c(all_sequences, Protein_sequence_fetch)
  }
  
  # Write all fetched sequences to the specified fasta file
  writeLines(all_sequences, fasta_file)
  
  return(paste("Protein Sequences written into", fasta_file))
}

# Example of usage:
# Using gene UIDs gathered with GeneRetrieval
# The sequences can be written in any format in the final file (FASTA format is recommended)
# batch_size must be smaller than 50
# ProteinFetcher(protein$PROT_UID,"protein_seq.fasta", batch_size = 40)
# ProteinFetcher("5594","protein_seq.fasta", batch_size = 40)

#====================================================================================================================#
#====================================================================================================================#

EnsemblRetrieval <- function(biomart_name,
                             dataset_name,
                             entrez_gene_ids,
                             host = "https://fungi.ensembl.org",
                             retrieve_gene = TRUE,
                             retrieve_mRNA = TRUE,
                             retrieve_protein = TRUE) {
  
  # Helper function to filter invalid characters
  filter_valid_sequences <- function(sequences, valid_chars) {
    gsub(paste0("[^", valid_chars, "]"), "", sequences)
  }
  
  # Connect to the specified Ensembl mart
  ensembl_mart <- useMart(host = host, biomart = biomart_name)
  
  # Select the dataset
  dataset <- useDataset(dataset_name, mart = ensembl_mart)
  
  # Get the Ensembl gene IDs
  results <- getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id"),
    filters = "entrezgene_id",
    values = entrez_gene_ids,
    mart = dataset
  )
  
  if (nrow(results) == 0) {
    stop("No results found for the provided Entrez gene IDs.")
  }
  
  gene_tibble <- tibble()
  mRNA_tibble <- tibble()
  protein_tibble <- tibble()
  
  # Helper function to retrieve and process sequences
  retrieve_sequences <- function(attributes, output_file, column_order, description_field, sequence_type) {
    sequences <- getBM(
      attributes = attributes,
      filters = "ensembl_gene_id",
      values = results$ensembl_gene_id,
      mart = dataset
    )
    
    sequence_tibble <- as_tibble(sequences) %>%
      dplyr::select(all_of(column_order)) %>%
      arrange(ensembl_gene_id)
    
    # Ensure sequences are correctly formatted
    valid_sequences <- sequences[[description_field]]
    valid_sequences <- valid_sequences[!is.na(valid_sequences) & valid_sequences != ""]
    
    if (sequence_type == "DNA") {
      valid_sequences <- filter_valid_sequences(valid_sequences, "ACGTNacgtn")
      fasta_sequences <- DNAStringSet(setNames(valid_sequences,
                                               paste(sequences$ensembl_gene_id,
                                                     sep = "_")))
    } else if (sequence_type == "RNA") {
      valid_sequences <- filter_valid_sequences(valid_sequences, "ACGUNacgun")
      fasta_sequences <- RNAStringSet(setNames(valid_sequences,
                                               paste(sequences$ensembl_gene_id,
                                                     sep = "_")))
    } else if (sequence_type == "AA") {
      valid_sequences <- filter_valid_sequences(valid_sequences, "ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy*")
      fasta_sequences <- AAStringSet(setNames(valid_sequences,
                                              paste(sequences$ensembl_gene_id,
                                                    sep = "_")))
    } else {
      stop("Invalid sequence type specified. Must be 'DNA', 'RNA', or 'AA'.")
    }
    
    writeXStringSet(fasta_sequences, filepath = output_file)
    
    return(sequence_tibble)
  }
  
  # Retrieve gene sequences
  if (retrieve_gene) {
    gene_tibble <- retrieve_sequences(
      attributes = c("ensembl_gene_id",
                     "entrezgene_id",
                     "description",
                     "gene_exon_intron"),
      output_file = "ensembl_gene_sequence.fasta",
      column_order = c("ensembl_gene_id", "entrezgene_id", "description"),
      description_field = "gene_exon_intron",
      sequence_type = "DNA"
    )
    cat("Gene sequences have been successfully retrieved and saved to FASTA files.\n")
  }
  
  # Retrieve mRNA sequences
  if (retrieve_mRNA) {
    mRNA_tibble <- retrieve_sequences(
      attributes = c("ensembl_gene_id",
                     "entrezgene_id",
                     "ensembl_transcript_id",
                     "description",
                     "transcript_exon_intron"),
      output_file = "ensembl_mRNA_sequence.fasta",
      column_order = c("ensembl_gene_id", "entrezgene_id", "ensembl_transcript_id", "description"),
      description_field = "transcript_exon_intron",
      sequence_type = "RNA"
    )
    cat("mRNA sequences have been successfully retrieved and saved to FASTA files.\n")
  }
  
  # Retrieve protein sequences
  if (retrieve_protein) {
    protein_tibble <- retrieve_sequences(
      attributes = c("ensembl_gene_id",
                     "ensembl_peptide_id",
                     "description",
                     "peptide"),
      output_file = "ensembl_protein_sequence.fasta",
      column_order = c("ensembl_gene_id", "ensembl_peptide_id", "description"),
      description_field = "peptide",
      sequence_type = "AA"
    )
    cat("Protein sequences have been successfully retrieved and saved to FASTA files.\n")
  }
  
  if (!retrieve_gene && !retrieve_mRNA && !retrieve_protein) {
    stop("No sequence type specified. Please select at least one sequence type to retrieve.")
  }
  
  return(list(gene_sequences = gene_tibble, mRNA_sequences = mRNA_tibble, protein_sequences = protein_tibble))
}

#====================================================================================================================#
#====================================================================================================================#

EnsemblInfo <- function(kingdom = "vertebrates") {
  # Define the base URLs for different kingdoms
  ensembl_urls <- list(
    vertebrates = "https://www.ensembl.org",
    plants = "https://plants.ensembl.org",
    fungi = "https://fungi.ensembl.org",
    protists = "https://protists.ensembl.org",
    metazoa = "https://metazoa.ensembl.org"
  )
  
  # Check if the specified kingdom is valid
  if (!kingdom %in% names(ensembl_urls)) {
    stop("Invalid kingdom specified. Valid options are: vertebrates, plants, fungi, protists, metazoa.")
  }
  
  # Get the BioMart list for the specified kingdom
  ensembl_url <- ensembl_urls[[kingdom]]
  bioMarts <- listMarts(host = ensembl_url)
  
  # Function to get datasets for a chosen BioMart
  getDatasetsForBioMart <- function(bioMartName) {
    selectedMart <- useMart(biomart = bioMartName, host = ensembl_url)
    datasets <- listDatasets(selectedMart)
    return(datasets)
  }
  
  # Create a tibble to store the results
  results <- tibble()
  
  # Loop through each BioMart and get its datasets
  for (i in 1:nrow(bioMarts)) {
    bioMartName <- bioMarts[i, "biomart"]
    datasets <- getDatasetsForBioMart(bioMartName)
    
    # Create a tibble for the current BioMart
    bioMartTibble <- tibble(
      biomart = bioMartName,
      dataset = datasets$dataset,
      host = ensembl_url,
      description = datasets$description
    )
    
    # Append to the results tibble
    results <- bind_rows(results, bioMartTibble)
  }
  
  return(results)
}


