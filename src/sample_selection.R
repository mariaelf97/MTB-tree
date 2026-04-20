


# sample from ncbi nucleotide samples, 10 per sub-lineage
ncbi_lineage <- fread("mnt/tb_seqs/nucleotide_db_isolates_lineages.tsv") 
IDs1 <- ncbi_lineage %>% select(Isolate,coll2014,freschi2020)%>%
  distinct(Isolate, .keep_all = TRUE)


#high quality assemblies
old_lineage <- fread("mnt/tb_seqs/all_isolates_lineages.tsv")
IDs2 <- old_lineage %>% select(accession_id,coll2014_lineage,freschi202_lineage) %>%
  distinct(accession_id, .keep_all = TRUE)

all_IDs<- rbind(IDs1,IDs2, use.names=FALSE)

set.seed(10)
ncbi_sample <- all_IDs %>% group_by(coll2014) %>% slice_sample(n=10) %>%
  distinct(Isolate, .keep_all = TRUE)
ncbi_sample <- ncbi_sample %>%
  mutate(lineage = str_extract(coll2014, "^[^.]+"))

ncbi_sample %>% write_tsv("mnt/tb_seqs/selected_sequences.tsv")
ncbi_sample %>% ungroup() %>% select(Isolate) %>% 
  write_tsv("mnt/tb_seqs/selected_sequences.txt")
