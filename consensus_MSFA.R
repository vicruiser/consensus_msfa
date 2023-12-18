## load necessary packages
requiredPackages = c('ggplot2', 'ggnewscale', 'cowplot',
                     'seqinr','stringr', 'data.table',
                     "dplyr", "tidyr", "Biostrings",
                     'grid','gridExtra','GetoptLong',
                     'ggseqlogo') 

suppressMessages(
  for (p in requiredPackages) {
    # if (!require(p, character.only = TRUE))
    #    install.packages(p)
    library(p, character.only = TRUE)
  }
)

#### functions
# map global msa to local msa (processed by FrustraEvo)
map_positions_msa = function(seq_ali, seq_subali) {
  # remove gaps from input sequences
  s1 = c2s(str_remove(seq_ali, '-'))
  s2 = c2s(str_remove(seq_subali, '-'))
  # check the position where the subalignment sequences starts
  loc = str_locate(s1, s2)
  # create a data frame containing AAs and the position of the protein sequence in the
  # global MSA
  sub_df = data.frame(aa = s2c(s2),
                      prot_pos = seq(loc[, 'start'],
                                     length(s2c(s2)) + loc[, 'start'] - 1,
                                     1))
  sub_df$sub_prot_pos = 1:nrow(sub_df)
  # create a data frame containing AAs and the position of the protein sequence in the
  # subMSA
  submsa_pos_df = data.frame (AA = seq_subali,
                              subali_pos = 1:length(seq_subali))
  
  submsa_pos_df = cbind(subset(submsa_pos_df, AA != '-'), sub_df)
  
  
  ddf = data.frame(aa = s2c(s1),
                   prot_pos = 1:length(s2c(s1)))
  
  msa_pos_df = data.frame (aa = seq_ali,
                           ali_pos = 1:length(seq_ali))
  
  msa_pos_df = cbind(subset(msa_pos_df, aa != '-'), ddf)
  
  
  mapping_df = right_join(submsa_pos_df, msa_pos_df[, c(1, 2, 4)], by = c("prot_pos", "aa"))
  
  return(mapping_df)
}


map_pdb_to_ali = function(pdb_seq, prot_seq, frust_res, seqid) {
  loc = str_locate(c2s(pdb_seq), c2s(prot_seq))
  if (any(is.na(loc))) {
    PWA  = pairwiseAlignment(subject = c2s(pdb_seq),
                             pattern = c2s(prot_seq))
    loc[, 'start'] = PWA@subject@range@start
    loc[, 'end'] = PWA@subject@range@width
  }
  ff = frust_res[loc[, "start"]:loc[, "end"],]
  colnames(ff)[colnames(ff) == "Res"] = "pdb_pos"
  ff$sub_prot_pos = 1:nrow(ff)
  ff$seqid = seqid
  return(ff)
}


mapping_all = function(msas, msa_glob, frust_files, input_dirs) {
  # list to store mapping data frames
  frustic_seqic_list = list()
  seqic = list()
  ## iterate over msas to get the SRFI indexes
  for (m in 1:length(msas)) {
    msa = msas[[m]]
    # get seqic
    slog = ggseqlogo(unlist(getSequence(getFrag(msa, 1,length(msa[[1]])), as.string = T)),
                     ncol = 100,
                     font= "akrobat_bold")+
           theme(axis.text.x = element_blank(),
            axis.text.y = element_text())+
            scale_x_continuous(breaks = c(1,seq(5,2998, by = 5)), expand = c(0, 0))
    
    LOGO_DATA = ggplot_build(slog)$data[[1]]
    LOGO_DATA$sub_prot_pos = round(LOGO_DATA$x)
    seqic[[m]] = LOGO_DATA %>% group_by(sub_prot_pos) %>% summarise(Entropy = max(y, na.rm = T))
    
    
    # for each sequence in the subalignment, get the mapping of positions
    # between the global and the local alignment and the frustratometeR results
    for (seqid in names(msa)) {
      # extract the sequence from the msa and submsa
      subali_seq = msa[str_detect(names(msa), seqid)][[1]]
      ali_seq = msa_glob[str_detect(names(msa_glob), seqid)][[1]]
      # get the mapping between msa and submsa for the seqid
      mm = map_positions_msa(ali_seq, subali_seq)
      mm$seqid = seqid
      mm$ref = seqrefs[m]
      # remove gaps from submsa sequence
      prot_seq = str_remove(subali_seq, '-')
      # read frustratometer results linked to seqid
      sel_frust_file = Sys.glob(file.path(input_dirs[m], 'Data', paste(seqid, '.done', sep=""),'FrustrationData', '*singleresidue'))
      frust_res = fread(sel_frust_file)
      # get the pdb sequence
      pdb_seq = frust_res$AA
      # mapping between submsa and frustratometer results
      ff = map_pdb_to_ali(pdb_seq,  prot_seq, frust_res, seqid)
      ff$ref = seqrefs[m]
      # join results
      m2 = left_join(mm, ff, by = c("sub_prot_pos", "ref", "AA", "seqid"))
      m2$cluster = cluster_names[m]
      # add mapping data frame to list of results
      frustic_seqic_list[[paste(m,seqid, sep ="_")]] = m2
    }
  }
  return(list(rbindlist(frustic_seqic_list), seqic))
}

# input arguments

GetoptLong(
  "msa_global=s@", "List of FrustraEvo output folders.",
  "input=s@", "Directory containing list of FrustraEvo output folders.",
  "cluster_names=s@", "Cluster names of each group. By default 'group_1', 'group_2',..., 'group_n'.",
  "width_plot=s@", "Width of the plot. Default = 80."
)
#
# global alignment
#msa_global = "/home/victoria/Escritorio/frustraheatmap/msas/Both.fasta"
msa_glob = read.fasta(msa_global, forceDNAtolower = F, seqtype = "AA")
#
# 2) frustraevo folders 
#input = Sys.glob("/home/victoria/Escritorio/frustraheatmap/example_data/FrustraEvo_*/pdb/FrustraEvo*") 
# 3) group names (optional)
#cluster_names = c("Alpha", "Beta")
if (length(cluster_names) < 1){
  cluster_names = paste('group', seq(1:length(input)), sep = "_")
}

# 4) max number of positions per row (default 100)
#width_plot = 80
## FALTA poner if no hay ná, default = 100

# Get paths to necesary files
input_dirs = Sys.glob(file.path(input,'*'))
msas_files = Sys.glob(file.path(input,'*', "OutPutFiles","MSA*.fasta"))
frustic_files = Sys.glob(file.path(input,'*', 'OutPutFiles', 'IC_SingleRes_*'))
frustic_files = frustic_files[!str_detect(frustic_files, '.txt')]
#seqic_files = Sys.glob(file.path(input,'*', 'OutPutFiles', 'SeqIC_*.tab'))
frust_files = Sys.glob(file.path(input,'*', 'Data','*.done','FrustrationData', '*.pdb_singleresidue'))

# load files
msas = lapply(msas_files, read.fasta, forceDNAtolower = F, seqtype = "AA")
frustic = lapply(frustic_files, fread)
seqic = list()#lapply(seqic_files, fread)

# get sequence of reference per protein group
seqrefs = unlist(lapply(frustic, function(x) unique(x$Prot_Ref)))

# get dataframe with positions of reference in msa and submsa
list_mapping = mapping_all(msas, msa_glob, frust_files, input_dirs)
df_mapping = list_mapping[[1]]
seqic = list_mapping[[2]]
# calculate median SRFI
median_srfi = df_mapping %>% group_by(cluster, ref, ali_pos, subali_pos) %>% summarise(median_SRFI = median(FrstIndex, na.rm=T))
# extract data for reference proteins
refs_srfi = subset(df_mapping, seqid %in% seqrefs ) 
# get consensus sequence in data frame
msas2 = lapply(msas_files, read.alignment, forceToLower = F, seqtype = "AA", format = "FASTA")
cmsa = lapply(msas2, consensus)
con = lapply(cmsa, function(c) data.frame(consensus=c, subali_pos = 1:length(c)) )
all <- do.call("rbind", con)
all$ref <- rep(seqrefs, lapply(con, nrow))
all$cluster <- rep(cluster_names, lapply(con, nrow))
# add consensus sequence to median SRFI calculations
con_median_srfi = left_join(median_srfi, all, by = c("subali_pos", "ref", "cluster"))
con_median_srfi = left_join(con_median_srfi, refs_srfi, by = c("ref", "ali_pos", "subali_pos", "cluster"))

# processing frustic data
frustic_df <- do.call("rbind", frustic)
frustic_df$cluster <- rep(cluster_names, lapply(frustic, nrow))
colnames(frustic_df)[colnames(frustic_df)=="Res"] = "sub_prot_pos"
colnames(frustic_df)[colnames(frustic_df)=="Num_Ref"] = "pdb_pos"
colnames(frustic_df)[colnames(frustic_df)=="Prot_Ref"] = "ref"


# processing of seqic data 
seqic_df <- do.call("rbind", seqic)
seqic_df$ref <- rep(seqrefs, sapply(seqic, nrow))
seqic_df$cluster <- rep(cluster_names, lapply(seqic, nrow))

#colnames(seqic_df)[colnames(seqic_df)=="Position"] = "sub_prot_pos"

# merge frustic and seqic data
frustic_seqic_df = merge(frustic_df, seqic_df, by = c("sub_prot_pos", "ref","cluster"))

# merge median SRFI, frustic and seqic data that contains the positions in the alignemnt
# needed to generate the consensus MSFA
plot_df = left_join(con_median_srfi, frustic_seqic_df, by = c("ref", "sub_prot_pos", "pdb_pos", "cluster"))
# determine the width per row of the plot
i = seq(0,max(plot_df$ali_pos), as.numeric(width_plot))+1
j = c(seq(as.numeric(width_plot) ,max(plot_df$ali_pos), as.numeric(width_plot)), max(plot_df$ali_pos))


# create the plot(s)
list_plots = list()
for (k in 1:length(i)){
  # subset data referring to just a certain width
  s =subset(plot_df, ali_pos %in% i[k]:j[k] & !is.na(Entropy))
  # generate plot
  g = ggplot(s,
             aes(x = ali_pos, y = cluster)) +
    geom_tile(
      aes(fill = median_SRFI),
      linejoin = "round",
      width = 0.96,
      height = 0.96,
      color = 'grey20'
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_text(aes(label  = consensus, size = Entropy), fontface = "bold") +
    coord_cartesian(clip = "off") +
    scale_fill_gradientn(
      values = c(1, .79, .5, 0.11, 0.0909, 0),
      colours = c("green", "grey", "grey", "grey", "red", "red"),
      name = "SRLF median index",
      limits = c(-1.2, 1.1),
      oob = scales::squish,
      na.value = "white"
    ) +
    scale_color_gradient2(
      mid = "grey",
      high = "black",
      midpoint = 0.5,
      limits = c(0, log2(20)),
      oob = scales::squish
    ) +
    theme_bw() +
    ylab(NULL) +
    xlab(NULL) +
    theme(legend.position = "top")
  # store plot
  list_plots[[k]] = g + theme(legend.position = "none")
}

# combine plots
combined_plot<-plot_grid(plotlist =list_plots,ncol=1)
# extract legend to have a single legend
legend1 <- get_legend(
  g +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "top")
)

# add legend to combined plot and legend using plot_grid()
p = plot_grid(NULL, legend1,NULL,   combined_plot, ncol=1,rel_heights = c(.1,.2,.1, 1))
# add x and y labels common for all the plots and not by row of the plot
y.grob <- textGrob("Cluster", 
                   gp=gpar(fontface="bold", col="grey20", fontsize=13), rot=90)
x.grob <- textGrob("Global MSA position", 
                   gp=gpar(fontface="bold", col="grey20", fontsize=13))
#pp = grid.arrange(arrangeGrob(p, left = y.grob, bottom = x.grob))

# save figure adjusting for height and width proportionally
height = length(cluster_names) + length(i) -1
svg("consensus_MSFA.svg", width = 15, height = height )
pp = grid.arrange(arrangeGrob(p, left = y.grob, bottom = x.grob))

dev.off()



# png('frustraheatmap.png', 
#     res=400,
#     height = 2000,
#     width = 6000)
# grid.arrange(arrangeGrob(p, left = y.grob, bottom = x.grob))
# dev.off()
# 
# 
# 
# fwrite(frustic_seqic, "frustic_seqic.txt")
