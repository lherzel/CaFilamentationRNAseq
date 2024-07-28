#!/usr/bin/env Rscript 

args = commandArgs(trailingOnly=TRUE)
readlength <- args[1]

if (length(readlength)==0) {
  readlength <- 101
}

# reads in Mio
counts1 <- as.numeric(read.table("counts1", header = F, stringsAsFactors = F)) / 1000000
counts2 <- as.numeric(read.table("counts2", header = F, stringsAsFactors = F)) / 1000000

classified <- read.table("classified.reads", header = F, stringsAsFactors = F)
names(classified) <- c("chr", "start.read", "end.read", "id.read", "block.count.read", "block.length.read", "block.start.read", 
                       "start.intron", "end.intron", "id.intron", "strand", "overlap", "category")

classified$length.intron <- classified$end.intron - classified$start.intron

print("input read")

types <- split(classified$category, classified$id.intron)
#types <- split(classified$type[1:1000000], classified$id.intron[1:1000000])

print("counting types of reads...")
types <- lapply(types, function(i) {
	i <- data.frame(table(i))
#	i$Var1 <- factor(i$Var1, levels = c("spliced", "unspliced", "internal"))
	return(i)
})

for (entry in 1:length(types)) {
  types[[entry]]$id.intron <- names(types)[entry]
}

types <- do.call("rbind", types)
row.names(types) <- NULL

print("adding intron length...")
classifiedRed <- unique(classified[,c("id.intron", "length.intron")])
types$length.intron <- classifiedRed$length.intron[match(types$id.intron, classifiedRed$id.intron, nomatch = NA)]

head(types)

print("normalizing internal read counts by length and library size..")
types$norm.byLength <- types$Freq / ((types$length.intron / (1000)) * counts1)

print("split by intron id to calculate intron retention fraction...")
types <- split(types, types$id.intron)
types <- lapply(types, function(i) {
  
  if (length(i$norm.byLength[i$i=="intronic"])>0) {
    internal <- i$norm.byLength[i$i=="intronic"]
    internal.abs <- i$Freq[i$i=="intronic"]
  } else {
    internal <- 0
    internal.abs <- 0
  }
  
  if (length(i$Freq[i$i=="unspliced"])>0) {
    unspliced <- i$Freq[i$i=="unspliced"]
  } else {
    unspliced <- 0
  }
  
  if (length(i$Freq[i$i=="spliced"])>0) {
    spliced <- i$Freq[i$i=="spliced"]
  } else {
    spliced <- 0
  }
  
  if (length(i$Freq[i$i=="alternative"])>0) {
    alternative <- i$Freq[i$i=="alternative"]
  } else {
    alternative <- 0
  }
  
	frcIR <- unspliced / (unspliced + 2*spliced)
	frcAS <- alternative / (alternative + spliced)
	frcIRwAS <- unspliced / (unspliced + 2*spliced + 2*alternative)
	frcASwUn <- alternative / (unspliced*0.5 + spliced + alternative)
	
	i <- data.frame(
	  id.intron = i$id.intron[1], 
	frcIR = frcIR,
	counts.frcIR = unspliced + spliced,
	frcAS = frcAS,
	counts.frcAS = alternative + spliced,
	frcIRwAS = frcIRwAS,
	frcASwUn = frcASwUn,
	count.spliced = spliced,
	count.unspliced = unspliced,
	count.intronic = internal.abs,
	count.intronic.rel = internal,
	count.alternative = alternative,
	counts.fwd = counts1,
	counts.rev = counts2,
	length.intron = i$length.intron[1])
	return(i)
})

types <- do.call("rbind", types)
anno <- unique(classified[,c("chr", "start.intron", "end.intron", "strand", "id.intron")])
anno$int <- paste(anno$chr, ":", anno$start.intron, "-", anno$end.intron, sep = "")

types <- merge(types, anno[,c("id.intron", "int", "strand")], by.x = "id.intron", by.y = "id.intron", all.x = T)

print("save table and exit")
write.table(types, "quantified", col.names = T, row.names = F, sep = "\t", quote = F)
