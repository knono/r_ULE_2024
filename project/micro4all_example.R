
#Load micro4all
library(micro4all)

##Load libraries
library(dada2); packageVersion("dada2")
library(ShortRead);packageVersion("ShortRead")
library(ggplot2);packageVersion("ggplot2")
library(tidyverse);packageVersion("tidyverse")



##SET PATH, FILES AND SAMPLE.NAMES
path <- "./rawData/" # CHANGE ME to the directory containing the fastq files
list.files(path) #list files in this folder

## SORT FORWARD AND REVERSE READS SEPARETELY; forward and reverse fastq filenames
## have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq 
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

##See name of the first sample (fnFs)
basename(fnFs[1])

##Apply strsplit
strsplit(basename(fnFs[1]), "_")

##COUNT NUMBER OF READS IN EACH SAMPLE BEFORE FILTERING 
raw_reads_count <- NULL #Creates an empty object

for (i in 1:length(fnFs)){ 
  #loop over fnFs to count number of sequences with length over readFastq
  raw_reads_count <- rbind(raw_reads_count,
                           length(ShortRead::readFastq(fnFs[i])))}

#Format table to give it sample.names
rownames(raw_reads_count)<- sample.names
colnames(raw_reads_count)<- "Number of reads"

#Get min and max number of sequences
min(raw_reads_count)

max(raw_reads_count)

##Histogram of sequences length

reads <- ShortRead::readFastq(fnFs) #store sequences in an R object

#Get lengths with unique
#get length accessing to width attribute within reads
uniques <- unique(reads@quality@quality@ranges@width) 

#Count number of sequences for each length
counts <- NULL #Creates a null object for storing counts
for (i in 1:length(uniques)) {
  counts<- rbind(counts,
                 length(which(reads@quality@quality@ranges@width==uniques[i])))
  
}

#format histogram table
histogram <-  cbind(uniques,counts)
colnames(histogram) <- c("Seq.length", "counts")

#Check histogram matrix
#Most sequences fall in expected sequence length
head(histogram[order(histogram[,1],decreasing = TRUE),])

# PLOT HISTOGRAM
hist(reads@quality@quality@ranges@width,
     main="Forward length distribution",
     xlab="Sequence length",
     ylab="Raw reads")

## VIEW AND SAVE QUALITY PLOT FOR FW AND RV ##
plotQualityProfile(fnFs[1:2])

plotQualityProfile(fnRs[1:2])

##FIGARO

##Creates a new folder called 'figaro' inside our path to store reads 
# after trimming for fígaro
figFs <- file.path(path, "figaro", basename(fnFs))
figRs <- file.path(path, "figaro", basename(fnRs))

##TRIMMING AT 295 pb
#inputFreads,outputFreads,inputRreads,outputRreads
out.figaro <- filterAndTrim(fnFs, figFs, fnRs, figRs, 
                            compress=TRUE, 
                            multithread=TRUE,
                            truncLen=c(295,295)) 

##RUN FIGARO

# figaro
# path to figaro program -i path to input - path to output -a length of your 
# amplicon without primers, -f primer forward length, -r primer reverse length

figaro <- system(("python3 /opt/figaro_py3.10/figaro/figaro.py \\
                  -i ./rawData/figaro \\
                  -o ./rawData/figaro \\
                  -a 426 -f 17 -r 21"),intern=TRUE)

head(figaro)

### FILTER AND TRIM SEQUENCES ACCORDING TO FIGARO ####
## Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(279,205),
                     maxN=0, maxEE=c(3,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, minLen=50)


head(out)

#### IDENTIFY PRIMERS ####

FWD <- "CCTACGGGNBGCASCAG"  ## CHANGE ME to your forward primer sequence
REV <- "GACTACNVGGGTATCTAATCC"  ## CHANGE ME to your reverse primer sequence

## VERIFY PRESENCE AND ORENTATION OF PRIMERS ##
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works with DNAString objects
                            # rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna),
               Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

## COUNT THE APPEARENCE AND ORIENTATION OF PRIMERS ##
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtFs[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = filtRs[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = filtFs[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = filtRs[[1]]))

## RUN CUTADAPT
cutadapt <- "/usr/bin/cutadapt" #path to cutadapt 

system2(cutadapt, args = c("--version")) # Run shell commands from R

##Create path to cutadapt sequences
path.cut <- file.path(path, "cutadapt") 

if(!dir.exists(path.cut)) dir.create(path.cut)

fnFs.cut <- file.path(path.cut, basename(filtFs))
fnRs.cut <- file.path(path.cut, basename(filtRs))


##Produce arguments for cutadapt. rc function creates the reverse complementary 
# of the provided sequence (FWD).  
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste0("-a", " ", "^",FWD,"...", REV.RC) 

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste0("-A"," ","^", REV, "...", FWD.RC)


# Run Cutadapt

for(i in seq_along(fnFs)) { system2 (cutadapt,
                                     args = c(R1.flags,
                                              R2.flags,
                                              "-n",2,
                                              "-m", 1,
                                              "--discard-untrimmed",
                                              "-j",0, 
                                              "-o", fnFs.cut[i],
                                              "-p", fnRs.cut[i], 
                                              filtFs[i], filtRs[i],
                                              "--report=minimal")) }  
  # -n 2 required to remove FWD and REV from reads
  # -m 1 is required to remove empty sequences for plotting quality plots
  # -j 0 automatically detect number of cores
  # -p output files
  # filtFs[i], filtRs[i], input files
  # Report minimal reports a summary

    
#Check primer presence after cutadapt
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

#### LEARN ERROR RATES ####
errF <- learnErrors(fnFs.cut, multithread=T, verbose=1) #

errR <- learnErrors(fnRs.cut, multithread=T, verbose=1) #

#Plot errors
plotErrors(errF, nominalQ=TRUE)

plotErrors(errR, nominalQ=TRUE)

#Sample inference
dadaFs <- dada(fnFs.cut, err=errF, multithread=TRUE)

dadaRs <- dada(fnRs.cut, err=errR, multithread=TRUE)

dadaFs[[1]]

#Set sample names
names(dadaFs) <- sample.names
names(dadaRs) <- sample.names

#Merge sequences
mergers <- mergePairs(dadaFs, fnFs.cut, dadaRs, fnRs.cut, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) 

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    method="consensus",
                                    multithread=TRUE,
                                    verbose=TRUE)

dim(seqtab.nochim)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.nochim)))  #Number of ASV of each length

#Get number of sequences per length
reads.per.seqlen <- tapply(colSums(seqtab.nochim),
                           factor(nchar(getSequences(seqtab.nochim))),
                           sum) #number of sequences for each length
reads.per.seqlen

## Plot length distribution
table_reads_seqlen <- data.frame(length=as.numeric(names(reads.per.seqlen)),
                                 count=reads.per.seqlen) #Create a table

ggplot(data=table_reads_seqlen, aes(x=length, y=count)) + geom_col() #plot it!

#Filter ASV length
seqtab.nochim <- seqtab.nochim[,
                               nchar(colnames(seqtab.nochim)) %in% seq(402,428)]

#Check number of sequences at each step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input",
                     "filtered",
                     "denoisedF",
                     "denoisedR",
                     "merged",
                     "nonchim")
rownames(track) <- sample.names


#Let's apply a mutate to get percentage, making it easier to analyze.
track <- track %>%as.data.frame() %>% mutate(Perc.filtered=filtered*100/input,
                                             Perc.denoisedF= denoisedF*100/filtered,
                                             Perc.denoisedR= denoisedR*100/filtered,
                                             Perc.merged= merged*100/filtered,
                                             Perc.nonchim= nonchim*100/merged,
                                             Perc.retained=nonchim*100/input
)

head(track)

write.table(track, file="dada2_track.txt", sep="\t")

#Classification with RDP
taxa_rdp <- assignTaxonomy(seqtab.nochim, "~/rawData/rdp_train_set_18_H.fa",
                           multithread=TRUE)

#Combine taxonomy and ASV abundance
ASV <- seqtab.nochim
ASVt <- t(ASV)

## SUBSTITUTE 'NA' WITH 'UNCLASSIFIED'and remove species column
taxa_rdp_na <- apply(taxa_rdp,2, tidyr::replace_na, "unclassified")[,-7]

##Create names for each ASV: if there are 1000 ASVs, call them 
# from ASV0001 to ASV1000
number.digit<- nchar(as.integer(nrow(ASVt)))
names <- paste0("ASV%0", number.digit, "d") #As many 0 as digits
ASV_names<- sprintf(names, 1:nrow(ASVt))

## CREATE AND SAVE ASV TABLE
ASV_table_classified<- cbind(as.data.frame(taxa_rdp_na,
                                           stringsAsFactors = FALSE),
                             as.data.frame(ASV_names, stringsAsFactors = FALSE),
                             as.data.frame(ASVt,stringsAsFactors = FALSE))

##Make rownames a new column, in order to keep sequences during the filtering 
# process
ASV_seqs <- rownames(ASV_table_classified)
rownames(ASV_table_classified) <- NULL
ASV_table_classified <- cbind(ASV_seqs, ASV_table_classified)

write.table(ASV_table_classified, file="ASV_table_classified.txt", sep="\t")

#MOCK community cleaning.
ASV_filtered_MOCK<-MockCommunity(ASV_table_classified,
                                 mock_composition,
                                 ASV_column = "ASV_names")
                                # mock_composition is a data frame included as 
                                # data in this package. See documentation for 
                                # details.

###CLEAN according to Bokulich####
# Get number of sequences of each ASV 
ASV_sums <- rowSums(ASV_table_classified[,9:ncol(ASV_table_classified)])
# Get total number of sequences
sum.total<-sum(ASV_sums)
# Apply percentage to sequence number
nseq_cutoff<-(0.005/100)*sum.total
# Filter table.
ASV_filtered_bokulich<- ASV_table_classified[which(ASV_sums>nseq_cutoff),]
# Sort table in ascending order of ASV names
ASV_table_bokulich<-ASV_filtered_bokulich[order(ASV_filtered_bokulich[["ASV_names"]]),]

#### CLEAN CHLOROPLAST SEQUENCES ####
ASV_final<-ASV_filtered_MOCK[(which(ASV_filtered_MOCK$Genus!="Streptophyta"
                                    & ASV_filtered_MOCK$Genus!="Chlorophyta"
                                    &ASV_filtered_MOCK$Genus!="Bacillariophyta"
                                    &ASV_filtered_MOCK$Family!="Streptophyta"
                                    & ASV_filtered_MOCK$Family!="Chlorophyta"
                                    &ASV_filtered_MOCK$Family!="Bacillariophyta"
                                    & ASV_filtered_MOCK$Family!="Mitochondria"
                                    & ASV_filtered_MOCK$Class!="Chloroplast"
                                    & ASV_filtered_MOCK$Order!="Chloroplast"
                                    & ASV_filtered_MOCK$Kingdom!="Eukaryota"
                                    & ASV_filtered_MOCK$Kingdom!="unclassified")),]

#Check Cyanobacteria
ASV_final[which(ASV_final$Phylum=="Cyanobacteria/Chloroplast"),]  


#Remove cyanobacteria at phylum level
ASV_final<-ASV_final[(which(ASV_final$Phylum!="Cyanobacteria/Chloroplast")),]

##Save table
write.table(ASV_final, file="ASV_final.txt", sep="\t")

#Analyze unclassified phyla
#Load libraries
library(seqinr); packageVersion("seqinr")
library(annotate); packageVersion("annotate")

#Get sequences of unclassified phyla
unclassified_phyla <-ASV_final[which(ASV_final$Phylum=="unclassified"),]
seq_unclassified_phyla<- as.list(unclassified_phyla$ASV_seqs)
#write fasta
write.fasta(seq_unclassified_phyla, names=unclassified_phyla$ASV_names,
            file.out="unclassified_phyla.fas",
            open = "w", nbchar =1000 , as.string = FALSE)

#Local BLASTn

# Para descargarlo de https://ftp.ncbi.nlm.nih.gov/blast/db/
# Actualmente 146 ficheros, por eso el contador antes del
# parallel
#
#> mkdir NCBI_nt_DB
#> cd NCBI_nt_DB
#> seq -w 00 146 | parallel wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt.{}.tar.gz
#
# Mas tarde descomprimirlo, ya vienen con formato; no hace falta indexarlo
#
# 1 #/bin/bash
# 2 for file in *.gz
# 3 do
# 4 tar -zxvpf "$file"
# 5 rm "$file"
# 6 done

# Una vez hecho, en mi caso
# > cd ~/Documentos/micro4all/project
# > ln -s /mnt/Datos/db/NBCI_nt_DB nt
#
# Para evitar sorpresas no uséis rm para borrar el enlace, en su caso
# usar unlink
# > unlink nt
#
# En mi caso voy a usar una copia de la base de datos cuando tenía 78 tar.gz
# ahora va por 146, así que paciencia.

# al final hice la consulta online, el fichero es pequeño

system("blastn -version")
system(("blastn -query unclassified_phyla.fas \\
        -remote -db nt \\
        -out unclassified_phyla_hits.txt \\
        -outfmt '6 std stitle' -show_gis \\
        -max_target_seqs 20 -parse_deflines"),
       intern = TRUE)
# -num_threads 10

####GET PHYLOGENETIC TREE####

library(phyloseq); packageVersion("phyloseq")
## GET THE FASTA FILE ##
ASVseqs_fasta<- as.list(ASV_final$ASV_seqs)
write.fasta(ASVseqs_fasta, names=ASV_final$ASV_names,
            file.out="ASV_tree.fas",
            open = "w",
            nbchar =1000 ,
            as.string = FALSE)

##Get the alignment ##
mafft <- "/usr/bin/mafft"     #path to program

system2(mafft, args="--version")
system2(mafft, args=c("--auto", "ASV_tree.fas>", "alignment"))

## Get the tree ##
FastTreeMP <- "/usr/bin/fasttreeMP"

system2(FastTreeMP, args="--version" )
system2(FastTreeMP, args = c("-gamma",
                             "-nt",
                             "-gtr",
                             "-spr",4 ,
                             "-mlacc", 2,
                             "-slownni",
                             "<alignment>",
                             "phy_tree"))# Run FastTreeMP

##Metadata from micro4all
metadata.micro <-micro4all::metadata
metadata.micro$samples <- str_replace_all(metadata.micro$samples,"-", "_")
rownames(metadata.micro) <- str_replace_all(rownames(metadata.micro),"-", "_")

####CREATE PHYLOSEQ OBJECT FROM TABLES

tax <- ASV_final[,2:8] #Tax
OTU <-  ASV_final[,9:ncol(ASV_final)] #ASV
colnames(OTU) <- str_replace_all(colnames(OTU),"-", "_")

#Sequences (BioString is imported with phyloseq)
dna<-Biostrings::DNAStringSet(ASV_final$ASV_seqs) 
names(dna)<- ASV_final$ASV_names

##ADD ASV NAMES
row.names(tax)<-ASV_final$ASV_names
row.names(OTU)<-ASV_final$ASV_names

#Check rownames are equal
identical(rownames(OTU), rownames(tax))


##Introduce phylogenetic tree
phy_tree <-  phyloseq::read_tree("phy_tree")
unrooted_tree<- phy_tree
ape::is.rooted(unrooted_tree)


##Produce root (we need a root for distance calculation)
tree_root<-ape::root(unrooted_tree, 1, resolve.root = T)
tree_root

ape::is.rooted(tree_root)

##CONVERT TO PHYLOSEQ FORMART
phy_OTUtable<-otu_table(OTU, taxa_are_rows = T)
phy_taxonomy<-tax_table(as.matrix(tax))
#Change here metadata.micro with mt (your own metadata object)
phy_metadata<-sample_data(metadata.micro) 

#Put everything into a phyloseq object
#Remove Tree_root when working with ITS
loc_phyloseq<-phyloseq(phy_OTUtable,phy_taxonomy,phy_metadata,dna,tree_root)

##Load packages
library(phyloseq); packageVersion("phyloseq")

library(Biostrings); packageVersion("Biostrings")

library(GUniFrac); packageVersion("GUniFrac")

library(phangorn); packageVersion("phangorn")

library(vegan); packageVersion("vegan")

library(pheatmap); packageVersion("pheatmap")

library(colorspace); packageVersion("colorspace")

#### RARECURVE ####
ASV <- as.data.frame(t(otu_table(loc_phyloseq)))
sample_names <- rownames(ASV)

#Generate rarefaction curves
rarecurve <- rarecurve(ASV, step = 100, label = F,)

#For each rarefaction curve, transform rarecurve output to a dataframe.
rarecurve_table <- lapply(rarecurve, function(x){
  b <- as.data.frame(x)
  b <- data.frame(ASV = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

#Produce a column with sample names and put everything in one data frame
names(rarecurve_table) <- sample_names
#for each rarecurve element, we produce a dataframe and bind them all together.
rarecurve_table <- purrr::map_dfr(rarecurve_table, function(x){ 
  z <- data.frame(x)
  return(z)
}, .id = "Sample")

#To color lines according to group, let's create a new column
#Coloring
color <- NULL
for (i in (1:length(rarecurve_table$Sample))) {
  ##Change metadata.micro to user mt and "Location" to variable of interest
  color <- rbind(color,
                 metadata.micro$location[which(metadata.micro$samples==rarecurve_table$Sample[i])])
}

#Bind this column
rarecurve_table_color <- cbind(rarecurve_table, color)
#Change Location to variable of interest
colnames(rarecurve_table_color) <- c(colnames(rarecurve_table),
                                     "Location")

## RARECURVE WITH GGPLOT ##

### Change location for your variable of interest and IMPORTANT TO GROUP BY SAMPLE
rareggplot<-ggplot(rarecurve_table_color,
                   aes(x=raw.read,
                       y=ASV,
                       colour=Location,
                       group=Sample)) + 
  theme_bw()+
  geom_point(aes(colour=Location),
             size=0.01)+ #Change location to your variable of interest
  ##Change location to your variable of interest, Change size for line thickness
  geom_line(aes(colour=Location),size=0.5)+ 
  geom_vline(aes(xintercept = min(sample_sums(loc_phyloseq))),
             lty=1,
             colour="black")+
  #Change location1, location2 and location3 by the values of your variable of 
  #interest (for example, control and inoculated)
  scale_fill_manual(values = c("location1"= "#33CC00",
                               "location2"= "#ffd624",
                               "location3"="#80CC88"))+ 
  #Change location1, location2 and location3 by the values of your variable of 
  #interest (for example, control and inoculated)
  scale_color_manual(values = c("location1"= "#33CC00",
                                "location2"= "#ffd624",
                                "location3"="blue"),
                     name="Location", #change Location to variable of interest
                     #change to values of your variable of interest
                     breaks=c("location1", "location2","location3"),
                     #change to values of your variable of interest
                     labels=c("Location 1","Location 2","Location 3"))+
  labs(title="Rarefaction curves", x="Number of sequences", y="Number of ASV")+
  guides(alpha="none")+
  theme(legend.key=element_blank(),
        legend.title.align = 0.85,
        legend.title = element_text(face="bold",size=14),
        axis.text = element_text(size=14),
        axis.title = element_text(size = 16),
        plot.title = element_text(hjust=0.5, face="bold", size=16),
        legend.text = element_text(size = 12))

rareggplot


#tiff
ggsave(filename = "rarefaction_curve.tiff",
       plot = rareggplot,device = tiff(),
       width = 18, height = 16, units = "cm", dpi = 800)

#svg
library(svglite)
ggsave(filename = "rarefaction_curve.svg", plot = rareggplot,device = svg())

#eps
postscript("curve_wito_RE.1.10.eps")
rareggplot
dev.off()  # Close the device

#Rarefy phyloseq
set.seed(10403)
rarefaction <- rarefy_even_depth(loc_phyloseq,
                                 sample.size = min(sample_sums(loc_phyloseq)),
                                 rngseed = FALSE)

#Estimate alpha indices and save it
alpha_diversity_rarefied <- estimate_richness(rarefaction,
                                              measures=c("InvSimpson",
                                                         "Shannon",
                                                         "Observed"))

## CALCULATE EVENNESS
Evenness_index <- as.data.frame(alpha_diversity_rarefied$Shannon/log(alpha_diversity_rarefied$Observed))
Evenness <- cbind(Evenness_index, rownames(alpha_diversity_rarefied))
colnames(Evenness) <- c("Evenness", "Samples")

#Generate final table
alpha_table <- cbind(alpha_diversity_rarefied,
                     Evenness$Evenness,
                     metadata.micro) #change metadata to user provided metadata
colnames(alpha_table)<- c("Observed",
                          "Shannon",
                          "InvSimpson",
                          "Evenness",
                          colnames(metadata.micro))#change metadata to user provided metadata

#Write final table on local machine
write.table(alpha_table, file="alpha_table.txt", sep="\t")

#### ALPHA PUB TABLE
#Select numeric columns and compute mean and SD
#Change location to variable of interest
alpha_mean <- aggregate(alpha_table[,1:4],
                        list(grouping=alpha_table$location),
                        mean)%>%  mutate_if(is.numeric,
                                            round, digits=2) 
alpha_sd <- aggregate(alpha_table[,1:4],
                      list(grouping=alpha_table$location),
                      sd)%>% mutate_if(is.numeric,
                                       round,
                                       digits=2)#Change location to variable of interest

#Paste mean ± SD 
mean_sd <- NULL
for (i in 2:5){
  mean_sd <- cbind(mean_sd,paste0(alpha_mean[,i], " +/- ", alpha_sd[,i]))}


#generate table and give it columns names
alpha_pub_table <- cbind(alpha_mean$grouping, mean_sd)
colnames(alpha_pub_table) <- c("Location",
                               "Observed richness",
                               "Shannon",
                               "InvSimpson",
                               "Evenness")#Change location to variable of interest

#Write table on your local machine
write.table(alpha_pub_table, file="alpha_pub_table.txt", sep="\t")

#### Apply levene and shapiro test ####

#change location to your variable of interest
levene<- levene.test.alpha(alpha_table, 4, "location") 

levene

#Write table on local machine
write.table(levene, file="levene_test.txt", sep="\t")

shapiro <- Shapiro(alpha_table, 4, "location") #change location to your variable of interest
shapiro

#Write table on local machine
write.table(shapiro, file="shapiro_test.txt", sep="\t")

#### Balanced anova ####
balanced_anova<- BalancedAnova(alpha_table,
                               numberOfIndexes = 4,
                               formula = "location") #change location to your variable of interest

balanced_anova[[1]]

#Write table on local machine
write.table(balanced_anova[[1]], file="balancedAnova_test.txt", sep="\t")

#### TUKEY ####
tukey <- Tukey.test(alpha_table, 4,
                    "location",
                    balanced=TRUE,)#change location to your variable of interest

tukey

#Write table on local machine
write.table(tukey, file="tukey_test.txt", sep="\t")

#### PRODUCE A TABLE FOR ALPHA PLOT ####
alpha_plot_table <- tidyr::pivot_longer(data = alpha_table,
                                        names_to = "Measure",
                                        values_to = "Value",
                                        cols=c(Observed,
                                               Shannon,
                                               InvSimpson,
                                               Evenness))

### ALPHA GRAPHIC GGPLOT ###
alpha_graphic <- ggplot(data = alpha_plot_table,
                        aes(x = location,
                            y = Value)) + #Change location to variable of interest
  facet_wrap(~factor(Measure,
                     levels=c("Observed",
                              "InvSimpson",
                              "Shannon",
                              "Evenness")),
             scale = "free") + #change levels order to reorder panels
  geom_boxplot()+
  theme(axis.text.x = element_text(size=13),
        legend.position="bottom",
        strip.text = element_text(size = 20),
        axis.text.y = element_text(size=15),
        axis.title.y=element_text(size=20)) +
  scale_x_discrete(limits=c("location1",
                            "location2",
                            "location3"),
                   breaks=c("location1",
                            "location2",
                            "location3"),   ##With breaks and labels you can change the name displayed on labels, #Change locations to values of variable of interest
                   labels=c("Location 1",#Change locations to values of variable of interest
                            "Location 2",
                            "Location 3")) +
  aes(fill=location)+
  scale_fill_manual(values = c( "location1"= "#33CC00",
                                #Change locations to values of variable of interest
                                "location2"= "#ffd624",
                                "location3"="blue"),
                    na.translate=FALSE) +
  theme(legend.key.size = unit(1, 'cm')) +
  ylab("Alpha Diversity Measure")

alpha_graphic

#use ggsave for saving ggplots
ggsave("alpha_plot.tiff",
       plot = alpha_graphic,
       width = 17,
       height = 30,
       units = "cm",
       dpi = 800)

library(edgeR); packageVersion("edgeR")

#### NORMALIZE ACCORDING TO EDGER ####
edgeR <- DGEList(counts = OTU,
                 samples = metadata.micro,
                 genes = tax)#Change metadata.micro for your own mt object
edgeR <- calcNormFactors(edgeR)

##EXTRACT NORMALIZED COUNTS
ASV_norm <- cpm(edgeR, normalized.lib.sizes=T, log=F)


##CREATE NORMALIZED PHYLOSEQ OBJECT
phy_OTU_norm<-otu_table(as.data.frame(ASV_norm,row.names=F), taxa_are_rows = T)
phy_taxonomy_norm<-tax_table(as.matrix(tax))
phy_metadata_norm<-sample_data(metadata.micro) #Change metadata.micro to your own mt file

##Add taxa names
taxa_names(phy_OTU_norm)<- taxa_names(phy_taxonomy_norm)
#Check
identical(rownames(ASV_norm), rownames(tax))
#> [1] TRUE

##Merge
norm_phyloseq<-phyloseq(phy_OTU_norm,
                        phy_taxonomy_norm,
                        phy_metadata_norm,
                        tree_root) #Remove tree_root when working with ITS

#### PERMANOVA ####
permanova<- Permanova(norm_phyloseq,
                      distances = c("bray",
                                    "unifrac",
                                    "wunifrac"),
                      formula = "location") #Change location to variable of interest

##Permanova result
permanova[[1]]

#Write table to local machine
write.table(permanova[[1]], file="permanova.txt",sep="\t")

#### BETADISPERSION ####
betadisper<- Betadispersion(loc_phyloseq,
                            distances = c("bray",
                                          "unifrac",
                                          "wunifrac"),
                            formula = "location") #change location to variable of interest

betadisper[1] #data frame with results for every distance method

#write file to your local machine
write.table(betadisper[1], file="betadisper.txt", sep="\t")


#### PAIRWISE ADONIS ####
pairwise<- PairwiseAdonisFun(norm_phyloseq,
                             distances = c("bray",
                                           "unifrac",
                                           "wunifrac"),
                             formula = "location",
                             pval=0.05) #Change location to your variable of interest


pairwise #data frame with all significant results

#write file to your local machine
write.table(pairwise, file="pairwise.txt", sep="\t")

#### BETA DIVERSITY PLOTS ####
var_list <-  c("bray", "unifrac", "wunifrac")
plot_type <-  c("PCoA", "NMDS")
combination_plot <- purrr::cross2(plot_type,var_list )


# Make plots.
plot_list = list()

for (i in 1:length(combination_plot)){
  pcoa <- ordinate(norm_phyloseq,combination_plot[[i]][[1]],combination_plot[[i]][[2]]) #get distance matrix
  plot <-  plot_ordination(norm_phyloseq, pcoa, type="samples", color="location")+ #change location to variable of interest (e.g., health status)
    geom_text(aes(label=location), hjust=0, vjust=0, show.legend=FALSE)+#change location variable of interest (e.g., health status)
    geom_point(size=4)
  if (combination_plot[[i]][[1]]=="NMDS"){ ## graphic details for NMDS
    plot <-  plot + xlab(paste("Axis 1"))+
      ylab(paste("Axis 2"))+
      theme(legend.position="bottom", plot.title = element_text(face="bold", hjust = 0.5), legend.title = element_blank(), legend.key = element_blank(),
            panel.border = element_rect(linetype = "solid", fill = NA),
            panel.background = element_rect(fill="white", colour="white"),
            panel.grid.major = element_line(colour="aliceblue", size=0.4),
            panel.grid.minor= element_line(colour = "aliceblue", size =0.4))+
      
      scale_color_manual(values=c("location1"= "#33CC00",
                                  "location2"= "#ffd624", "location3"="blue"))+#change locations to values of variable of interest (eg. healthy, infected)
      guides(color=guide_legend(nrow=2,byrow=TRUE))+
      guides(shape=guide_legend(nrow=2,byrow=TRUE))+
      theme(plot.title = element_text(face="bold", hjust = 0.5))+
      ggtitle(paste("Bacterial Community", combination_plot[[i]][[1]], "on", combination_plot[[i]][[2]], "distances", round(pcoa$stress, digits = 3)))
    
    plot_list[[i]] = plot
    
  }
  else {## graphic details for PCoA
    plot= plot +   xlab(paste("PCo 1", paste("(",round(pcoa$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
      ylab(paste("PCo 2", paste("(",round(pcoa$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
      
      theme(legend.position="bottom", plot.title = element_text(face="bold", hjust = 0.5), legend.title = element_blank(), legend.key = element_blank(),
            panel.border = element_rect(linetype = "solid", fill = NA),
            panel.background = element_rect(fill="white", colour="white"),
            panel.grid.major = element_line(colour="aliceblue", size=0.4),
            panel.grid.minor= element_line(colour = "aliceblue", size =0.4))+
      
      scale_color_manual(values=c("location1"= "#33CC00",
                                  "location2"= "#ffd624", "location3"="blue"))+#change locations to values of variable of interest (eg. healthy, infected)
      guides(color=guide_legend(nrow=2,byrow=TRUE))+
      guides(shape=guide_legend(nrow=2,byrow=TRUE))+
      theme(plot.title = element_text(face="bold", hjust = 0.5))+
      ggtitle(paste("Bacterial Community", combination_plot[[i]][[1]], "on", combination_plot[[i]][[2]], "distances"))
    
    plot_list[[i]] = plot
    
  }
  
}


#Visualize plot
plot_list[[1]]

#Save plot on local machine
ggsave("bray-curtis_pcoa.tiff", plot = plot_list[[1]],
       width = 17, height = 30, units = "cm", dpi = 800 )



####PCOA####
#Calculate distances
pcoa <- ordinate(norm_phyloseq,"PCoA","bray") #get distance matrix, change bray to choosen distance 

#Produce plot
pcoa_bray <-  plot_ordination(norm_phyloseq, pcoa, type="samples", color="location")+ #change location to variable of interest
  geom_text(aes(label=location), hjust=0, vjust=0, show.legend=FALSE)+#change location to variable of interest
  geom_point(size=4) + 
  xlab(paste("PCo 1", paste("(",round(pcoa$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
  ylab(paste("PCo 2", paste("(",round(pcoa$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
  theme(legend.position="bottom", plot.title = element_text(face="bold", hjust = 0.5), legend.title = element_blank(), legend.key = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        panel.background = element_rect(fill="white", colour="white"),
        panel.grid.major = element_line(colour="aliceblue", size=0.4),
        panel.grid.minor= element_line(colour = "aliceblue", size =0.4))+
  scale_color_manual(values=c("location1"= "#33CC00",
                              "location2"= "#ffd624", "location3"="blue"))+#change locations to values of variable of interest (eg. healthy, infected)
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  guides(shape=guide_legend(nrow=2,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("Bacterial Community PCoA on Bray Curtis distances") #change title according to distance used

#Visualize plot 
pcoa_bray

#Save it to local machine
ggsave("bray-curtis_pcoa.tiff", plot = pcoa_bray,width = 17,
       height = 30, units = "cm", dpi = 800 )

#### NMDS####
#Calculate distances
nmds <- ordinate(norm_phyloseq,
                 "NMDS",
                 "bray") #get distance matrix, change bray to choosen distance 

#Produce plot
nmds_bray <-  plot_ordination(norm_phyloseq, nmds, type="samples", color="location")+ #change location to variable of interest (e.g., health status)
  geom_text(aes(label=location), hjust=0, vjust=0, show.legend=FALSE)+#change location to variable of interest
  geom_point(size=4)+
  xlab(paste("Axis 1"))+
  ylab(paste("Axis 2"))+
  theme(legend.position="bottom", plot.title = element_text(face="bold", hjust = 0.5), legend.title = element_blank(), legend.key = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA),
        panel.background = element_rect(fill="white", colour="white"),
        panel.grid.major = element_line(colour="aliceblue", size=0.4),
        panel.grid.minor= element_line(colour = "aliceblue", size =0.4))+
  
  scale_color_manual(values=c("location1"= "#33CC00",
                              "location2"= "#ffd624", "location3"="blue"))+ #change locations to values of variable of interest (eg. healthy, infected)
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  guides(shape=guide_legend(nrow=2,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle(paste("Bacterial Community NMDS on Bray-Curtis distances", round(nmds$stress, digits = 3))) #change title according to distance used

#Visualize plot 
nmds_bray


#Save it to local machine
ggsave("bray-curtis_nmds.tiff",
       plot = nmds_bray,width = 17, height = 30, units = "cm", dpi = 800 )


#### ORDINATE PHYLOSEQ OBJECT ####
PCOA_bray  <-  ordinate(norm_phyloseq,
                        "PCoA",
                        "bray",
                        formula="location") #Change to your variable of interest
#Get number of axes
length(PCOA_bray$values$Relative_eig)
#> [1] 23
#Get plot
plot_scree(PCOA_bray)

##Get explained variance in first 10 axis
sum(PCOA_bray$values$Relative_eig[1:10]) #change 1:10 to select the number of axes

##CALCULATE RELAIVE ABUNDANCE, ASV, SAMPLES
sample_relabun_ASV <- transform_sample_counts(loc_phyloseq, function(x){x/sum(x)}*100)

##CONSTRUCT THE TABLE
OTU_sample <-  as.data.frame((otu_table(sample_relabun_ASV)))
taxonomy_sample <- as.data.frame(tax_table(sample_relabun_ASV))
identical(rownames(OTU_sample),rownames(taxonomy_sample))

sample_relabun_ASV_table  <- cbind(taxonomy_sample, OTU_sample)

sample_relabun_ASV_table <- sample_relabun_ASV_table[sort(sample_relabun_ASV_table$ASV_names,
                                                          decreasing = FALSE),]

#Check relative abundance sums 100
colSums(sample_relabun_ASV_table[,8:ncol(sample_relabun_ASV_table)])

##See table
head(sample_relabun_ASV_table)

#Save ASV sample table on your local machine
write.table(sample_relabun_ASV_table, file="ASV_sample_relabun.txt", sep="\t")

### GET TABLE BY SAMPLE AT GENUS LEVEL ###
#Glom phyloseq
loc_phyloseq_genus <- tax_glom(loc_phyloseq, taxrank = "Genus") #Change Genus to the desired taxonomic rank 
sample_relabun_genus <- transform_sample_counts(loc_phyloseq_genus, function(x){x/sum(x)}*100)
##Extract elements from phyloseq
OTU_sample_genus <-  as.data.frame((otu_table(sample_relabun_genus)))
taxonomy_sample_genus  <- as.data.frame(tax_table(sample_relabun_genus)[,-7])
identical(rownames(OTU_sample_genus),rownames(taxonomy_sample_genus))

sample_table_genus  <- cbind(taxonomy_sample_genus, OTU_sample_genus)

##Agglomerate unclassified in otu table

#Get sum of unclassified sequences
sample_table_genus_unclass <- sample_table_genus %>% subset(Genus=="unclassified", select=c(7:ncol(sample_table_genus))) %>%
  colSums() %>%  t() %>% as.data.frame() %>% cbind(Kingdom="unclassified", Phylum="unclassified", Class="unclassified",
                                                   Order="unclassified", Family="unclassified", Genus="unclassified",  .)

#Bind sample_table_genus wito unclassified with sample_table_genus_class
sample_table_genus_final <- rbind(subset(sample_table_genus, Genus!="unclassified"), sample_table_genus_unclass)

#Check relative abundance sums 100
colSums(sample_table_genus_final[,8:ncol(sample_table_genus_final)])

#Sort table from most to least abundant genera (calculated as sum of abundance between samples)
sample_table_genus_final <- sample_table_genus_final[order(rowSums(sample_table_genus_final[,7:ncol(sample_table_genus_final)]), decreasing=TRUE),]

#See first lines of this table
head(sample_table_genus_final)

#Save to your local machine
write.table(sample_table_genus_final, file="sample_table_genus.txt", sep="\t")

### GROUP GENUS TABLE BY LOCATION

#Save genus abundances
otu_genus <- sample_table_genus_final[,7:ncol(sample_table_genus_final)] %>% t() %>% as.data.frame()
#Save taxonomy data
tax_genus <- sample_table_genus_final[,1:6]

#Calculate OTU mean abundance based on grouping factor (e.g., location)
location_mean_genus <- aggregate(otu_genus, by=list(metadata.micro$location), FUN=mean)%>% column_to_rownames("Group.1") %>% t() #Change metadata.micro and location to your metadata and variable of interest
#Calculate OTU SD  based on grouping factor (e.g., location) and change colnames
location_SD_genus <- aggregate(otu_genus, by=list(metadata.micro$location), FUN=sd)%>% column_to_rownames("Group.1")  %>% t()  %>% #Change metadata.micro and location to your metadata and variable of interest
  as.data.frame() %>% rename_with(.fn= ~paste0(colnames(location_mean_genus), "SD"))

#Merge mean abundance, SD and taxonomy.
genus_location_table <- merge(tax_genus, location_mean_genus, by=0) %>%column_to_rownames("Row.names") %>%
  merge(location_SD_genus, by=0) %>% column_to_rownames("Row.names")

#Check abundances sum 100
colSums(genus_location_table[,7:ncol(genus_location_table)])

#Sort table from most to least abundant genera (calculated as sum of abundance between samples)
genus_location_table <- genus_location_table[order(rowSums(genus_location_table[,7:ncol(genus_location_table)]), decreasing=TRUE),]

#View table
head(genus_location_table)

#Save table on your local machine
write.table(genus_location_table, file="genus_abund_location.txt", sep="\t")

#Save OTU data
otu_location_ASV <- sample_relabun_ASV_table[,8:ncol(sample_relabun_ASV_table)] %>% t() %>% as.data.frame()
#Save Taxnomy data
tax_location_ASV <- sample_relabun_ASV_table[,1:8]

#Calculate OTU mean abundance based on grouping factor (e.g., PLOT)
location_ASV_mean <- aggregate(otu_location_ASV, by=list(metadata.micro$location), FUN=mean)%>% column_to_rownames("Group.1") %>% t() #Change metadata.micro and location to your metadata and variable of interest
#Calculate OTU SD  based on grouping factor (e.g., PLOT) and change colnames
location_ASV_SD <- aggregate(otu_location_ASV, by=list(metadata.micro$location), FUN=sd)%>% column_to_rownames("Group.1")  %>% t()  %>% #Change metadata.micro and location to your metadata and variable of interest
  as.data.frame() %>% rename_with(.fn= ~paste0(colnames(location_ASV_mean), "SD"))

#Merge mean abundance, SD and taxonomy.
ASV_location_table <- merge(tax_location_ASV, location_ASV_mean, by=0) %>%column_to_rownames("Row.names") %>%
  merge(location_ASV_SD, by=0) %>% column_to_rownames("Row.names")

#Check abundances sum 100
colSums(ASV_location_table[,8:ncol(ASV_location_table)])

#Save to local machine
write.table(ASV_location_table, file="ASV_abund_location.txt", sep="\t")

###TAXONOMICAL PROFILE GENUS LEVEL###
#Format table. 
#Add a column with the mean abundance between conditions and arrange it un decreasing order. 
sample_table_genus_sort <- genus_location_table %>% mutate(mean=rowMeans(genus_location_table[,7:ncol(genus_location_table)])) %>% arrange(desc(mean))

#We have to use pivot_longer to get a table in the right format for ggplot
#We filter out unclassified
taxonomic_genus_location <-  sample_table_genus_sort[1:21,] %>% filter(Genus!="unclassified") %>% dplyr::select(6:9) %>% pivot_longer(!Genus, names_to="Samples", values_to="Abundance")##we select 2:21 because the unclassified sequences are the most abundant.

library(gdata);packageVersion("gdata")

#Get labels for location
location_label <- levels(as.factor(unique(taxonomic_genus_location$Samples)))

#Get unique names for genera
unique_genera <- unique(as.vector(taxonomic_genus_location$Genus)) #Change Genus to Phylum when needed

#REORDER FACTORS LEVELS IN DATA FRAME
taxonomic_genus_location$Genus=reorder.factor(taxonomic_genus_location$Genus,new.order=rev(unique_genera))#Change Genus to Phylum when needed

sorted_labels_genus<- as.data.frame(unique_genera)

##CREATE AN EXPRESSION FOR GGPLOT WITH ITALICS WHEN NEEDED.
sorted_labels_ggplot <- sapply(sorted_labels_genus$unique_genera,
                               function(x) parse(text = paste0("italic('",as.character(x), "')")))

library(randomcoloR);packageVersion("randomcoloR")
library(scales);packageVersion("scales")

#Produce random color and visualize them
#Get number of genera
length(unique_genera) #This will be the number of genera to color

colors_genus <-  randomColor(count=length(unique_genera), hue=c("random"))
show_col(colors_genus)

#Modify colors
colors_genus[5] <-  "#448800"

#Produce a palette of colours
c24 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "pink", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "blue","orange", "darkgrey"
)

#### BACTERIAL GENERA TAXONOMIC PROFILE ####

ggplot_title <- "Bacterial genera by location"

ggGenera=ggplot(taxonomic_genus_location, aes(x=Samples, y=Abundance, fill=Genus, order=Genus)) + #Change genus to phylum
  geom_bar(stat="identity", position="stack")+
  scale_fill_manual(values=c24,#Change c24 to colors_genus if you are using ramdom colors
                    breaks=unique_genera, #Include genera names
                    labels=sorted_labels_ggplot)+ #Include italic tags
  
  labs(y="Mean relative abundance (%)", x=NULL, fill="Genus",title=ggplot_title)+
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits=location_label,
                   labels=location_label)+
  scale_y_continuous(expand=c(0.01,0.01),
                     breaks=c(0,10,20,30,40,50,60,70),
                     labels=c("0","10", "20","30","40","50","60","70"), 
                     limits = c(NA, 70))+
  theme_bw()+
  theme(panel.border = element_rect(colour="black"), axis.title.x=element_blank(),
        plot.title = element_text(face="bold", hjust = 0.5, size=16), axis.text = element_text(size = 14),
        axis.text.x = element_text(face="bold", size=10, angle = 45, vjust=1, hjust = 1), axis.title.y = element_text(size = 16),
        legend.key.size = unit(0.9, "cm"), legend.text = element_text(size = 8),
        legend.title = element_text(size=7, face="bold"), legend.title.align=0.5)


#Visualize plot
ggGenera

#Save it
ggsave("genera_profile.tiff", plot = ggGenera,width = 35, height = 30, units = "cm", dpi = 800 )

library(microbiome);packageVersion("microbiome")

#### ANCOMBC COMPARISONS ####
ANCOM_location_genus <- ancomloop(input_object_phyloseq = loc_phyloseq, grouping = "location", ancom.options = list(global=FALSE, struc_zero=TRUE, n_cl=8),out.unclassified = TRUE, tax.level="Genus") #When performing ancombc at genus level, we can filter unclassified genera out with out.unclassified and tax.level arguments #REMEMBER! Change location1 according to your metadata

#Visualize first lines 
head(ANCOM_location_genus[["location1"]])

#Save any of the table on your local machine
write.table(ANCOM_location_genus[["location1"]], file="ancom_location1_genus.txt", sep="\t")

###Filter ANCOM Results
#Get results for location1
ANCOM_location1_genus <- ANCOM_location_genus[["location1"]] #Change location1 to your variable of interest
#Get only significant results
ANCOM_loc1vsloc2_genus_sig <- ANCOM_location1_genus[,1:19] %>%  filter(location1vslocation2_diff_abn,TRUE) #Change [,1:19] according to the comparison of interest

#Write results on your local machine
write.table(ANCOM_loc1vsloc2_genus_sig, file="ANCOM_loc1loc2_sig.txt", sep="\t")

#Bind phylum name to genus name
#Create a new variable for graphics
graphic_loc1vsloc2_genus<-ANCOM_loc1vsloc2_genus_sig

#Loop to add a new column called Classification and paste phylum and genus names
for (i in 1:nrow(graphic_loc1vsloc2_genus)){
  graphic_loc1vsloc2_genus$Classification[i]<- paste(graphic_loc1vsloc2_genus$Phylum[i],graphic_loc1vsloc2_genus$Genus[i],sep="|")
}

#Get the last classified taxonomy level
for (i in 1:nrow(graphic_loc1vsloc2_genus)){ #Loop over rows
  if (isTRUE(graphic_loc1vsloc2_genus$Genus[i]=="unclassified")& isTRUE(graphic_loc1vsloc2_genus$Family[i]!="unclassified")){
    graphic_loc1vsloc2_genus$Classification[i] <- paste0("Family", " ", graphic_loc1vsloc2_genus$Family[i],
                                                         "|", graphic_loc1vsloc2_genus$Genus[i])}
  else if (isTRUE(graphic_loc1vsloc2_genus$Family[i]=="unclassified") &isTRUE(graphic_loc1vsloc2_genus$Order[i]!="unclassified")){
    graphic_loc1vsloc2_genus$Classification[i] <- paste0("Order", " ", graphic_loc1vsloc2_genus$Order[i],
                                                         "|", graphic_loc1vsloc2_genus$Genus[i])
  }
  else if (isTRUE(graphic_loc1vsloc2_genus$Order[i]=="unclassified") &isTRUE(graphic_loc1vsloc2_genus$Class[i]!="unclassified")){
    graphic_loc1vsloc2_genus$Classification[i] <- paste0("Class", " ", graphic_loc1vsloc2_genus$Class[i],
                                                         "|", graphic_loc1vsloc2_genus$Genus[i])
  }
  else if (isTRUE(graphic_loc1vsloc2_genus$Class[i]=="unclassified") &isTRUE(graphic_loc1vsloc2_genus$Phylum[i]!="unclassified")){
    graphic_loc1vsloc2_genus$Classification[i] <- paste0("Phylum", " ", graphic_loc1vsloc2_genus$Phylum[i],
                                                         "|", graphic_loc1vsloc2_genus$Genus[i])
  }
  else if (isTRUE(graphic_loc1vsloc2_genus$Phylum[i]=="unclassified")){
    graphic_loc1vsloc2_genus$Classification[i] <- paste0("Phylum", " ", "unclassified")
  }
  else if (isTRUE(graphic_loc1vsloc2_genus$Genus[i]!="unclassified")) {
    graphic_loc1vsloc2_genus$Classification[i] <- paste0("Phylum", " ", graphic_loc1vsloc2_genus$Phylum[i],
                                                         "|", graphic_loc1vsloc2_genus$Genus[i])
    
  }}

#Filter out unclassified phyla
graphic_loc1vsloc2_genus <- graphic_loc1vsloc2_genus %>% filter(Phylum!="unclassified") 

#Get the last classified taxonomy level for ANCOM at ASV level
for (i in 1:nrow(graphic_loc1vsloc2_ASV)){
  if (isTRUE(graphic_loc1vsloc2_ASV$Genus[i]=="unclassified")& isTRUE(graphic_loc1vsloc2_ASV$Family[i]!="unclassified")){
    graphic_loc1vsloc2_ASV$Classification[i] <- paste0("Family", " ", graphic_loc1vsloc2_ASV$Family[i],
                                                       "|", graphic_loc1vsloc2_ASV$Genus[i], "|", graphic_loc1vsloc2_ASV$ASV_names[i])}
  else if (isTRUE(graphic_loc1vsloc2_ASV$Family[i]=="unclassified") &isTRUE(graphic_loc1vsloc2_ASV$Order[i]!="unclassified")){
    graphic_loc1vsloc2_ASV$Classification[i] <- paste0("Order", " ", graphic_loc1vsloc2_ASV$Order[i],
                                                       "|", graphic_loc1vsloc2_ASV$Genus[i],"|", graphic_loc1vsloc2_ASV$ASV_names[i])
  }
  else if (isTRUE(graphic_loc1vsloc2_ASV$Order[i]=="unclassified") &isTRUE(graphic_loc1vsloc2_ASV$Class[i]!="unclassified")){
    graphic_loc1vsloc2_ASV$Classification[i] <- paste0("Class", " ", graphic_loc1vsloc2_ASV$Class[i],
                                                       "|", graphic_loc1vsloc2_ASV$Genus[i],"|", graphic_loc1vsloc2_ASV$ASV_names[i])
  }
  else if (isTRUE(graphic_loc1vsloc2_ASV$Class[i]=="unclassified") &isTRUE(graphic_loc1vsloc2_ASV$Phylum[i]!="unclassified")){
    graphic_loc1vsloc2_ASV$Classification[i] <- paste0("Phylum", " ", graphic_loc1vsloc2_ASV$Phylum[i],
                                                       "|", graphic_loc1vsloc2_ASV$Genus[i],"|", graphic_loc1vsloc2_ASV$ASV_names[i])
  }
  else if (isTRUE(graphic_loc1vsloc2_ASV$Phylum[i]=="unclassified")){
    graphic_loc1vsloc2_ASV$Classification[i] <- paste0("Phylum", " ", "unclassified")
  }
  else if (isTRUE(graphic_loc1vsloc2_ASV$Genus[i]!="unclassified")) {
    graphic_loc1vsloc2_ASV$Classification[i] <- paste0("Phylum", " ", graphic_loc1vsloc2_ASV$Phylum[i],
                                                       "|", graphic_loc1vsloc2_ASV$Genus[i],"|", graphic_loc1vsloc2_ASV$ASV_names[i])
    
  }}

##Create the data frame for representation, with Classification, Log fold change and standard deviation
logfold_df <-  graphic_loc1vsloc2_genus %>% dplyr::select(Classification, location1vslocation2_lfc, location1vslocation2_se)

colnames(logfold_df)<- c("Genus", "LogFold_location1vslocation2", "SD")

##Filter out 0 log fold change and assign name according to the direction of change
logfold_df  <-  logfold_df %>%
  filter(LogFold_location1vslocation2 != 0) %>% arrange(desc(LogFold_location1vslocation2)) %>%
  mutate(group = ifelse(LogFold_location1vslocation2 > 0, "Location2", "Location1")) #Change location to your variable of interest from your mt object

logfold_df$Genus = factor(logfold_df$Genus, levels = logfold_df$Genus)

write.table(file="waterfall_table.txt", logfold_df, sep="\t")

#FILTER WHEN NEEDED. This can be used to filter the results according to log fold change level
range(logfold_df$LogFold_location1vslocation2)

logfold_df_filtered <- rbind(logfold_df[which(logfold_df$LogFold_location1vslocation2>=1.5),], logfold_df[which(logfold_df$LogFold_location1vslocation2<=-1.5),])

##Waterfall plot
waterfall_location1 <-  ggplot(data = logfold_df_filtered,#If not filtering, use logfold_df in plotting
                               aes(x = Genus, y = LogFold_location1vslocation2, fill = group, color = group)) + #Change Genus to taxonomic rank of interest and logFold_location1vslocation2 to your object
  geom_bar(stat = "identity", width = 0.7,
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = LogFold_location1vslocation2 - SD, ymax = LogFold_location1vslocation2 + SD), width = 0.2, #Change LogFold_location1vslocation2 to your object
                position = position_dodge(0.05), color = "black") +
  labs(x = NULL, y = "Log fold change",
       title = "Waterfall Plot for the Location Effect") + #change title 
  theme_bw() +
  theme(legend.position = "right",
        legend.title =element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1), plot.margin=margin(10,10,10,100))

#Visualize plot
waterfall_location1

#Save it
ggsave("waterfall_location1.tiff", plot = waterfall_location1,width = 17, height = 30, units = "cm", dpi = 800 )

#Again, bind phylum name and genus name 
pyramidal_df<- genus_location_table
for (i in 1:nrow(pyramidal_df)){
  pyramidal_df$Classification[i]<- paste(pyramidal_df$Phylum[i],pyramidal_df$Genus[i],sep="|")
}
#Filter to include only significant genera according to ANCOM
pyramidal_loc1vsloc2_filt <- pyramidal_df[which(pyramidal_df$Classification %in%logfold_df_filtered$Genus ),]

#Sort it in the same way as waterfall plot
pyramidal_loc1vsloc2_filt <- pyramidal_loc1vsloc2_filt[match(logfold_df_filtered$Genus,pyramidal_loc1vsloc2_filt$Classification),]

#Now, transform it to ggplot format
pyramidal_loc1vsloc2 <- pyramidal_loc1vsloc2_filt %>% dplyr::select(c(13, 7:8)) %>% #select columns with Classification and abundance
  pivot_longer(!Classification, names_to="Samples", values_to="Abundance")#

write.table(file="pyramidal_table.txt", pyramidal_loc1vsloc2, sep="\t")

#Split table according to location
loc1pyrm <- pyramidal_loc1vsloc2[which(pyramidal_loc1vsloc2$Samples=="location1"),]#Change location1 to variable of interest
loc2pyrm <- pyramidal_loc1vsloc2[which(pyramidal_loc1vsloc2$Samples=="location2"),]#Change location1 to variable of interest

#Create vector with factors
loc1pyrm$Classification <- factor(loc1pyrm$Classification, levels=loc1pyrm$Classification)
loc2pyrm$Classification <- factor(loc2pyrm$Classification, levels=loc2pyrm$Classification)

#Get limits for plotting
min(pyramidal_loc1vsloc2$Abundance)

max(pyramidal_loc1vsloc2$Abundance)

#### PYRAMIDAL PLOT ####
pyramidalggplot=ggplot(data = pyramidal_loc1vsloc2, mapping= aes(x = Classification, y=Abundance, fill = Samples), colour="white")+
  geom_bar(loc1pyrm, stat="identity", mapping=aes(y=-Abundance))+ #in data, give the first split dataframe
  geom_bar(data=loc2pyrm, stat="identity")+ #in data, give the second split dataframe
  theme_bw()+
  scale_y_continuous(expand=c(0,0), labels=abs, limits=c(-4,4), breaks=seq(-4,4,1))+ #Change limits here
  labs(y="Relative abundance (%)")+
  theme(legend.title=element_blank(), axis.title.y=element_blank(), axis.text.y= element_text(size = 8, face="bold"),
        axis.text.x = element_text(size=14, face="bold"), axis.title.x = element_text(size=18, face="bold"), legend.key.size = unit(1.1, "cm"),legend.text = element_text(size = 16), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank())+
  coord_flip()


#Visualize plot
pyramidalggplot

#Save it
ggsave("pyramidal_loc1_loc2.tiff", plot = pyramidalggplot,width = 17, height = 30, units = "cm", dpi = 800 )

#Get genus location table in reverse order
genus_location_table_rev <- genus_location_table[nrow(genus_location_table):1,]
#Again, bind phylum name and genus name 
pyramidal_df_abun<- genus_location_table_rev
for (i in 1:nrow(pyramidal_df_abun)){
  pyramidal_df_abun$Classification[i]<- paste(pyramidal_df_abun$Phylum[i],pyramidal_df_abun$Genus[i],sep="|")
}
#Filter to include only significant genera according to ANCOM
pyramidal_loc1vsloc2_filt_abund <- pyramidal_df_abun[which(pyramidal_df_abun$Classification %in%logfold_df_filtered$Genus ),]

#Now, transform it to ggplot format
pyramidal_loc1vsloc2_abund <- pyramidal_loc1vsloc2_filt_abund %>% dplyr::select(c(13, 7:8)) %>% #select columns with Classification and abundance
  pivot_longer(!Classification, names_to="Samples", values_to="Abundance")#

#Split table according to location
loc1pyrm_abund <- pyramidal_loc1vsloc2_abund[which(pyramidal_loc1vsloc2_abund$Samples=="location1"),]#Change location1 to variable of interest
loc2pyrm_abund <- pyramidal_loc1vsloc2_abund[which(pyramidal_loc1vsloc2_abund$Samples=="location2"),]#Change location1 to variable of interest

#Create vector with factors
loc1pyrm_abund$Classification <- factor(loc1pyrm_abund$Classification, levels=loc1pyrm_abund$Classification)
loc2pyrm_abund$Classification <- factor(loc2pyrm_abund$Classification, levels=loc2pyrm_abund$Classification)

#### PYRAMIDAL PLOT ####
pyramidalggplot_abund=ggplot(data = pyramidal_loc1vsloc2_abund, mapping= aes(x = Classification, y=Abundance, fill = Samples), colour="white")+
  geom_bar(loc1pyrm_abund, stat="identity", mapping=aes(y=-Abundance))+ #in data, give the first split dataframe
  geom_bar(data=loc2pyrm_abund, stat="identity")+ #in data, give the second split dataframe
  theme_bw()+
  scale_y_continuous(expand=c(0,0), labels=abs, limits=c(-4,4), breaks=seq(-4,4,1))+ #Change limits here
  labs(y="Relative abundance (%)")+
  theme(legend.title=element_blank(), axis.title.y=element_blank(), axis.text.y= element_text(size = 8, face="bold"),
        axis.text.x = element_text(size=14, face="bold"), axis.title.x = element_text(size=18, face="bold"), legend.key.size = unit(1.1, "cm"),legend.text = element_text(size = 16), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank())+
  coord_flip()


#Visualize plot
pyramidalggplot_abund

#Save it
ggsave("pyramidal_loc1_loc2.tiff", plot = pyramidalggplot,width = 17, height = 30, units = "cm", dpi = 800 )

### PREPARE DATA FOR CAP #### 
#Read phyloseq and get relative abundance
cap_phyloseq <- readRDS(file = "cap_phyloseq")
cap_phyloseq_rel <- transform_sample_counts(cap_phyloseq,function(x){x/sum(x)}*100)

#Get ASV and metadata
ASV_cap<-as.data.frame(t(otu_table(cap_phyloseq_rel)))
metadata_cap <- as.matrix(sample_data(cap_phyloseq))
metadata_cap <- as.data.frame(metadata_cap)

#metadata_cap <- metadata_cap[1:ncol(metadata_cap)-1] #Remove column Season_NucleicAcid

#Compute CAP1 and CAP0 models
CAP1<-capscale(ASV_cap~Water+ln_P+SOM+N+ph_KCl+CN+Na_Exch+Clay+Sand+Slime+Texture, data=metadata_cap,distance = "bray")#Add here, after ASV_cap, all your parameters

CAP0<-capscale(ASV_cap~1, data=metadata_cap, distance = "bray")

#Get formula with ordistep
CAP_ordi<-ordistep(CAP0, scope=formula(CAP1))

#### ENVFIT ####
ef <-  envfit (CAP1~ph_KCl, data = metadata_cap, perm = 999) # AFter CAP1, introduce the significant parameters

#Now we can adjust pvalues of envfit
ef.adj <- ef
pvals.adj <- p.adjust (ef$vectors$pvals, method = 'BH')
ef.adj$vectors$pvals <- pvals.adj
ef.adj

#### CAP PLOT ####
#First, get distance matrix
distance_matrix <- phyloseq::distance(physeq = cap_phyloseq_rel, method = "wunifrac")

#Then, generate ordination with CAP (that is, constrained ordination to ph_KCl)
CAP_wunifrac <-  ordinate(physeq = cap_phyloseq_rel, method = "CAP",distance = distance_matrix,
                          formula = ~ph_KCl+Clay)#Change parameters to yours

#Produce plot
CAP_wunifrac_plot  <- plot_ordination(physeq = cap_phyloseq_rel, ordination = CAP_wunifrac,
                                      color = "Species", axes = c(1,2)) +
  aes(shape = Species) +
  geom_point(aes(colour = Species), alpha = 1, size = 3.5) +
  scale_shape_manual(values=c("Pinaster"=16,"Other"=15), #Change Pinaster and other according to your variable of interest
                     breaks=c("Pinaster","Other"))+#Change Pinaster and other according to your variable of interest
  scale_color_manual(values = c("Pinaster"="brown", "Other"="orange"), #Change Pinaster and other according to your variable of interest
                     name="Host genotype", 
                     breaks=c("Pinaster", "Other"),
                     labels=expression(paste("Confirmed ", italic("P. pinaster")), "Pine forest samples"))+ #This is used to set the legend, change it according to your preference
  guides(shape="none")+
  ggtitle("CAP on Weighted Unifrac distance")+
  theme_bw()+
  theme(legend.key=element_blank(),
        legend.title.align = 0.85,
        legend.title = element_text(face="bold"),
        axis.text = element_text(size=14),
        axis.title = element_text(size = 16),
        plot.title = element_text(hjust=0.5, face="bold"),
        legend.text = element_text(size = 16))+
  geom_hline(aes(yintercept = c(0.00)), lty=2, colour="grey")+
  geom_vline(aes(xintercept = c(0.00)), lty=2, colour="grey")

#Define arrows aesthetic and labels
arrowmat <-  vegan::scores(CAP_wunifrac, display = "bp")
arrowdf <-  data.frame(labels = rownames(arrowmat), arrowmat) 
arrow_map  <-  aes(xend = CAP1, yend = CAP2,x = 0, y = 0, shape = NULL, color = NULL, label=labels)
label_map <- aes(x = 1.1 * CAP1, y = 1.1 * CAP2, shape = NULL,color = NULL, label=labels)
arrowhead  <-  arrow(length = unit(0.02, "npc"))

#Introduce them to plot
CAP_wunifrac_plot +
  geom_segment(mapping = arrow_map, size = c(ph_KCl=1,Clay=1), data = arrowdf, color = 
                 c(ph_KCl="black",Clay="black"),arrow = arrowhead) +#Change pH and Clay to your physycochemical variables of interest
  geom_text(mapping = label_map, size = 4,data = arrowdf, show.legend = FALSE)

##Load library

library(car);packageVersion("car")
library(stats);packageVersion("stats")
library(BiocGenerics); packageVersion("BiocGenerics")
library(agricolae);packageVersion("agricolae")
library(pairwiseAdonis);packageVersion("pairwiseAdonis")
library(dunn.test);packageVersion("dunn.test")


#### ALPHA DIVERSITY WITHOUT MICRO4ALL ####
##Balanced anova
res.aov_observed <-  aov(Observed ~ location, data = alpha_table)#Change OBSERVED for index of interest and location to mt variable
summary(res.aov_observed)

##Unbalanced anova
res.aov_Observed_u <-  Anova(res.aov_observed, type = "III")# Change res.aov_observed for every index (res.aov_shannon, etc...) 
res.aov_Observed_u 

##POST-HOC TEST, TUKEYHSD 
post_observed <- TukeyHSD(res.aov_observed, which = "location")$location # Change res.aov_observed for every index (res.aov_shannon, etc...) and location to mt variable
which(post_observed[,"p adj"]<0.05) 

##Post-hoc, dunn Test
dunn.result<- dunn.test(x=alpha_table$Observed, g=alpha_table$location, method = "BH")

##Levene test and shapiro test
leveneTest(Observed ~ location, data = alpha_table)#Change OBSERVED for index of interest  and location to mt variable

#Shapiro
aov_residuals_observed = residuals(object = res.aov_observed) # Change res.aov_observed for every index (res.aov_shannon, etc...)
shapiro.test(x = aov_residuals_observed) 

##Kruskal Wallis
k_observed <- kruskal.test(Observed ~ location, data = alpha_table) #Change OBSERVED for index of interest and location to mt variable
k_observed

#Wilcoxon test
pairwise.wilcox.test(alpha_table$Observed, alpha_table$location, p.adjust.method = "BH") #Change OBSERVED for index of interest and location to mt variable

#### BETA DIVERSITY WITHOUT MICRO4ALL ####
#PERMANOVA with bray curtis distances
df <- data.frame(sample_data(norm_phyloseq))
bray <- phyloseq::distance(norm_phyloseq, method="bray") #Change bray to include another distance method
permanova_all <- adonis2(bray~ location, data=df)  #Change location to mt variable of interest
permanova_all

#Betadisper with bray curtis distances
disper_d <- betadisper(bray, df$location) #Change location to mt variable of interest
permutest(disper_d) # NOT SIGNIFICANT 

#Pairwise adonis
d <-  data.frame(sample_data(norm_phyloseq))
table_vegan <- cbind(t(OTU), df) 
ncolum <- ncol(t(OTU))

## PAIRWISE ADONIS ON BRAY CURTIS ## 
pw_result <- pairwise.adonis(table_vegan[,1:ncolum], table_vegan$location, p.adjust.m = "BH") #Change location to mt variable of interest
pw_result_sig <- pw_result[which(pw_result$p.adjusted<0.05),]
pw_result_sig # None

