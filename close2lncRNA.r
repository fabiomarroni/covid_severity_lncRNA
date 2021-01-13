# Run with --help flag for help.
# Modified 12/30/2018 by Fabio Marroni

suppressPackageStartupMessages({
  require(optparse)
})

option_list = list(
  make_option(c("-I", "--input_file"), type="character", default="COVID19_HGI_B1_ALL_20201020.txt.gz_1.0E-5.txt",
              help="Input file (cuffdiff output) [default= %default]", metavar="character"),
  make_option(c("-G", "--gff_file"), type="character", default="GCF_000001405.39_GRCh38.p13_genomic.gff",
              help="Input file (cuffdiff output) [default= %default]", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="ncRNA_COVID19_HGI_B1_ALL_20201020.txt.gz_1.0E-5.txt", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$input_file)) {
  stop("WARNING: No input file specified with '-I' flag.")
} else {  cat ("Input file is ", opt$input_file, "\n")
  input_file <- opt$input_file  
  }

if (is.null(opt$gff_file)) {
  stop("WARNING: No gff file specified with '-G' flag.")
} else {  cat ("gff file is ", opt$gff_file, "\n")
  gff_file <- opt$gff_file  
  }

  if (is.null(opt$out)) {
  stop("WARNING: No input directory specified with '-I' flag.")
} else {  cat ("Output file is ", opt$out, "\n")
  outfile <- opt$out  
  }


  
close2lncRNA<-function(input_file, gff_file, outfile) 
{
	library(data.table)
	#I implemented this on windows, so I couldn't exploit the power of data.table(cmd=grep...) and I used this workaround to give gff a decent form
	gff<-scan(gff_file,what="",sep="\n",comment.char = "#")
	#Throw away a lot of minor scaffolds: we don't want them now
	gff<-gff[grep("^NC_",gff)]
	#Trnasform the gff in a 9 column matrix as it SHOULD be!
	mgff<-matrix(unlist(strsplit(gff,"\t")),ncol=9,byrow=T)
	#gff has its onw chromosome names, we need just a quick trick to bring them to the usual name
	mgff[,1]<-gsub("NC_0000","",unlist(lapply(strsplit(mgff[,1],"\\."),"[",1)))
	mgff[,1][mgff[,1]=="NC_012920"]<-"MT"
	mgff[,1][mgff[,1]=="024"]<-"Y"
	mgff[,1][mgff[,1]=="023"]<-"X"
	mgff[,1]<-gsub("0","",mgff[,1])
	mgff<-data.frame(mgff,stringsAsFactors=F)
	#Convert start and end positions to numeric
	mgff[,4]<-as.numeric(mgff[,4])
	mgff[,5]<-as.numeric(mgff[,5])
	#To simplify our check, we only keep exon. We can retrieve any info from this
	mgff<-mgff[mgff[,3]=="exon",]
	#Now that we are fine with the gff we can read the result file
	gwares<-fread(input_file,data.table=F)
	gwares$lncRNA<-"No"
	for(aaa in 1:nrow(gwares))
	{
	sgff<-mgff[as.character(mgff[,1])==gwares$"#CHR"[aaa],]
	check_start<-sgff[which.min(abs(sgff[,4]-gwares$POS[aaa])),]
	check_end<-sgff[which.min(abs(sgff[,5]-gwares$POS[aaa])),]
	#If the closest exon (be it exon start or exon end) to the SNP is a ncRNA, write down the name
	#I guess I should be proud of this shit for retrieving the lncRNA name!
	if(length(grep("ncRNA",check_start))>0) gwares$lncRNA[aaa]<-unlist(lapply(strsplit(unlist(lapply(strsplit(check_start$X9,"gene="),"[",2)),",|;"),"[",1))
	if(length(grep("ncRNA",check_end))>0) gwares$lncRNA[aaa]<-unlist(lapply(strsplit(unlist(lapply(strsplit(check_end$X9,"gene="),"[",2)),",|;"),"[",1))
	}
	write.table(gwares,outfile,row.names=F,quote=F,sep="\t")	
	browser()
}
close2lncRNA(input_file=input_file,gff_file=gff_file,outfile=outfile)