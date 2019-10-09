library(maftools)
args=commandArgs(T)

transformSegments = function(segmentedData, build = "hg19"){
  build.opts = c('hg19', 'hg18', 'hg38')

  if(!build %in% build.opts){
    stop('Available reference builds: hg18, hg19, hg38')
  }

  if(build == 'hg19'){
    chr.lens = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                 155270560, 59373566)
  } else if(build == 'hg18'){
    chr.lens = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
                 158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
                 114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
                 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  } else if(build == 'hg38'){ #hg38
    chr.lens = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  } else{
    stop('Available reference builds: hg18, hg19, hg38')
  }

  segmentedData[,Start_Position := as.numeric(as.character(Start_Position))]
  segmentedData[,End_Position := as.numeric(as.character(End_Position))]

  #Replace chr x and y with numeric value (23 and 24) for better ordering
  segmentedData$Chromosome = gsub(pattern = 'chr', replacement = '', x = segmentedData$Chromosome, fixed = TRUE)
  segmentedData$Chromosome = gsub(pattern = 'X', replacement = '23', x = segmentedData$Chromosome, fixed = TRUE)
  segmentedData$Chromosome = gsub(pattern = 'Y', replacement = '24', x = segmentedData$Chromosome, fixed = TRUE)

  segmentedData$Chromosome = factor(x = segmentedData$Chromosome, levels = 1:24, labels = 1:24)

  segmentedData = segmentedData[order(Chromosome, Start_Position, decreasing = FALSE)]

  seg.spl = split(segmentedData, segmentedData$Chromosome)

  seg.spl.transformed = seg.spl[[1]]
  if(nrow(seg.spl.transformed) > 0){
    seg.spl.transformed$Start_Position_updated = seg.spl.transformed$Start_Position
    seg.spl.transformed$End_Position_updated = seg.spl.transformed$End_Position
  }
  chr.lens.sumsum = cumsum(chr.lens)
  
  for(i in 2:length(seg.spl)){
  
    x.seg = seg.spl[[i]]
    if(nrow(x.seg) > 0){
      x.seg$Start_Position_updated = x.seg$Start_Position + chr.lens.sumsum[i-1]
      x.seg$End_Position_updated = x.seg$End_Position + chr.lens.sumsum[i-1]
    }
    seg.spl.transformed = rbind(seg.spl.transformed, x.seg, fill = TRUE)
  }

  return(seg.spl.transformed)
}
setGeneric(name = "getCytobandSummary", function(x) standardGeneric("getCytobandSummary"))
setMethod(f = "getCytobandSummary",signature = "GISTIC", function(x) x@cytoband.summary)

rungiastcmain <- function (gistic = NULL, fdrCutOff = 0.1, markBands = NULL, color = NULL, 
    ref.build = "hg19", loh=NULL, cytobandOffset = 0.01, txtSize = 0.8, 
    cytobandTxtSize = 0.6) 
{
    g = getCytobandSummary(gistic)
    g = g[qvalues < fdrCutOff]
    g[, `:=`(Chromosome, sapply(strsplit(x = g$Wide_Peak_Limits, 
        split = ":"), "[", 1))]
    g[, `:=`(loc, sapply(strsplit(x = g$Wide_Peak_Limits, split = ":"), 
        "[", 2))]
    g[, `:=`(Start_Position, sapply(strsplit(x = g$loc, split = "-"), 
        "[", 1))]
    g[, `:=`(End_Position, sapply(strsplit(x = g$loc, split = "-"), 
        "[", 2))]
    g.lin = transformSegments(segmentedData = g[, .(Chromosome, 
        Start_Position, End_Position, qvalues, Cytoband, Variant_Classification)],build=ref.build)
    if (is.null(color)) {
        color = c(Amp = "red", Del = "blue")
    }
    gis.scores = transformSegments(segmentedData = gistic@gis.scores, 
        build = ref.build)
    gis.scores$amp = ifelse(test = gis.scores$Variant_Classification == 
        "Del", yes = -gis.scores$G_Score, no = gis.scores$G_Score)
    gis.scores$ystart = ifelse(test = gis.scores$Variant_Classification == 
        "Del", yes = -cytobandOffset, no = cytobandOffset)
    fdrCutOff = -log10(fdrCutOff)
    gis.scores$Variant_Classification = ifelse(test = as.numeric(gis.scores$fdr) > 
        fdrCutOff, yes = gis.scores$Variant_Classification, no = "neutral")
    gis.scores$Variant_Classification = factor(gis.scores$Variant_Classification, 
        levels = c("neutral", "Amp", "Del"))
    if (ref.build == "hg19") {
        chr.lens = c(249250621, 243199373, 198022430, 191154276, 
            180915260, 171115067, 159138663, 146364022, 141213431, 
            135534747, 135006516, 133851895, 115169878, 107349540, 
            102531392, 90354753, 81195210, 78077248, 59128983, 
            63025520, 48129895, 51304566, 155270560, 59373566)
    }
    else if (ref.build == "hg18") {
        chr.lens = c(247249719, 242951149, 199501827, 191273063, 
            180857866, 170899992, 158821424, 146274826, 140273252, 
            135374737, 134452384, 132349534, 114142980, 106368585, 
            100338915, 88827254, 78774742, 76117153, 63811651, 
            62435964, 46944323, 49691432, 154913754, 57772954)
    }
    else if (ref.build == "hg38") {
        chr.lens = c(248956422, 242193529, 198295559, 190214555, 
            181538259, 170805979, 159345973, 145138636, 138394717, 
            133797422, 135086622, 133275309, 114364328, 107043718, 
            101991189, 90338345, 83257441, 80373285, 58617616, 
            64444167, 46709983, 50818468, 156040895, 57227415)
    }
    else {
        stop("ref.build can only be hg18, hg19 or hg38")
    }
    chr.lens.cumsum = cumsum(chr.lens)
    nchrs = length(unique(gis.scores$Chromosome))
    chr.labels = c(1:22, "X", "Y")
    chr.tbl = data.table::data.table(chr = chr.labels, start = c(1, 
        chr.lens.cumsum[1:length(chr.lens.cumsum) - 1]), end = chr.lens.cumsum)
    chr.tbl$color = rep(c("black", "white"), length = nrow(chr.tbl))
    y_lims = pretty(gis.scores[, amp], na.rm = TRUE)
    gis.scores$Variant_Classification = factor(x = gis.scores$Variant_Classification, 
        levels = c("neutral", "Amp", "Del"))
    gis.scores = split(gis.scores, as.factor(as.character(gis.scores$Variant_Classification)))
    par(mar = c(2, 4, 2, 1))
    plot(NA, NA, xlim = c(0, chr.lens.cumsum[length(chr.lens.cumsum)]), 
        ylim = range(y_lims), axes = FALSE, xlab = NA, ylab = NA)
    abline(v = chr.tbl$end, h = y_lims, lty = 2, col = grDevices::adjustcolor("gray70", 
        0.25))
    axis(side = 2, at = round(y_lims, 2), las = 2)
    mtext(text = "G-Score", side = 2, line = 2.8, cex = 1.2)
    if (nrow(gis.scores[["neutral"]]) > 0) {
        segments(x0 = gis.scores[["neutral"]]$Start_Position_updated, 
            y0 = 0, x1 = gis.scores[["neutral"]]$End_Position_updated, 
            y1 = gis.scores[["neutral"]]$amp, col = "gray70")
    }
    if (nrow(gis.scores[["Amp"]]) > 0) {
        segments(x0 = gis.scores[["Amp"]]$Start_Position_updated, 
            y0 = 0, x1 = gis.scores[["Amp"]]$End_Position_updated, 
            y1 = gis.scores[["Amp"]]$amp, col = color[1])
    }
    if (nrow(gis.scores[["Del"]]) > 0) {
        segments(x0 = gis.scores[["Del"]]$Start_Position_updated, 
            y0 = 0, x1 = gis.scores[["Del"]]$End_Position_updated, 
            y1 = gis.scores[["Del"]]$amp, col = color[2])
    }
    rect(xleft = chr.tbl$start, ybottom = -cytobandOffset, xright = chr.tbl$end, 
        ytop = cytobandOffset, col = chr.tbl$color)
    text(y = 0, x = apply(chr.tbl[, 2:3], 1, mean), labels = chr.tbl$chr, 
        cex = cytobandTxtSize, col = c("white", "black"))
    gis.scores = data.table::rbindlist(l = gis.scores, use.names = TRUE, 
        fill = TRUE)
    if (is.null(markBands)) {
        markBands = g.lin[order(qvalues)][1:5, Cytoband]
    }
    if (all(length(markBands) == 1 & markBands == "all")) {
        mb = g.lin
    }
    else {
        mb = g.lin[Cytoband %in% markBands]
    }
    if (nrow(mb) == 0) {
        message("Available cytobands: ")
        print(getCytobandSummary(x = gistic)[qvalues < fdrCutOff])
        stop(paste("Could not find provided cytobands:", paste(markBands, 
            collapse = ", ")))
    }
    data.table::setkey(x = gis.scores, Chromosome, Start_Position_updated, 
        End_Position_updated)
    cyto_peaks_scores = data.table::foverlaps(y = gis.scores[, 
        .(Chromosome, Start_Position_updated, End_Position_updated, 
            amp)], x = mb[, .(Chromosome, Start_Position_updated, 
        End_Position_updated, Cytoband)])
    cyto_peaks_scores = cyto_peaks_scores[order(Cytoband, abs(amp), 
        decreasing = TRUE)][Cytoband %in% mb$Cytoband][!duplicated(Cytoband)]
    if (nrow(cyto_peaks_scores) > 1) {
        wordcloud::textplot(x = cyto_peaks_scores$Start_Position_updated, 
            y = cyto_peaks_scores$amp, words = cyto_peaks_scores$Cytoband, 
            new = FALSE, font = 3, cex = txtSize)
    }
    else {
        text(x = cyto_peaks_scores$Start_Position_updated, y = cyto_peaks_scores$amp, 
            labels = cyto_peaks_scores$Cytoband, font = 3, cex = txtSize)
    }
    if (!is.null(loh)) {
	for (p in loh$pos){
	    ypos=rnorm(1,mean=-0.7,sd=0.1)
	    xpos=rnorm(1,mean=p,sd=100000000)
	    text(x = xpos, y = ypos, labels = loh$gene[which(loh$pos == p)], col='red', font = 3, cex = txtSize)
	    lines(c(p,xpos),c(0,ypos-0.05),lty=2)
    	}
    }
}

laml.gistic = readGistic(gisticAllLesionsFile = args[1], gisticAmpGenesFile = args[2], gisticDelGenesFile = args[3], gisticScoresFile = args[4], isTCGA = FALSE)
sigloh=args[5]
ref_build=args[6]
outfile=args[7]

sigloh_data=read.table(sigloh,sep='\t',header=TRUE)
pdf(outfile,width=10,height=10)
if (nrow(sigloh_data) == 0){
	rungiastcmain(gistic = laml.gistic, markBands = "all",ref.build=ref_build)
}else{
	rungiastcmain(gistic = laml.gistic, markBands = "all",ref.build=ref_build,loh=sigloh_data)
}
dev.off()

