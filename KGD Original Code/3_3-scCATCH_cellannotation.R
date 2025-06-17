##scCATCH for cell annotation
# #install
# install.packages("scCATCH")
#scCATCH
library(scCATCH)

#load TN.combined
TN.combined <- readRDS(file = "output/Human_combined.rds")
#DefaultAssay(TN.combined) <- "RNA"

# # demo_geneinfo
# demo_geneinfo()
# revise gene symbols
data.input <- GetAssayData(TN.combined, assay = "RNA", slot = "data") # normalized data matrix
data.input <- rev_gene(data = data.input, data_type = "data", species = "Human", geneinfo = geneinfo)
#create scCATCH object with createscCATCH(). Users need to provide the normalized data and the cluster for each cell.
labels <- Idents(TN.combined)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
obj <- createscCATCH(data = data.input, cluster = as.character(meta$group))
# demo_geneinfo
demo_marker()
# The most strict condition to identify marker genes
obj <- findmarkergene(object = obj, 
                      species = "Human", 
                      marker = cellmatch, 
                      tissue = c('Adipose tissue','Blood','Peripheral blood','Bone','Cartilage','Subcutaneous adipose tissue',
                                 'Hair follicle','Lung','Muscle','Skin','Dermis','Lymph node','Lymphoid tissue',
                                 'Pluripotent stem cell','Skeletal muscle','Umbilical cord blood','Plasma',
                                 'Umbilical cord','Spleen','Serum','Bone marrow','Placenta','Embryonic stem cell','Kidney',
                                 'Pancreas','Pancreatic islet','Pyloric gland','Pancreatic acinar tissue'
                      ), 
                      use_method = "1")
obj <- findcelltype(object = obj)
obj@celltype
write.csv(obj@celltype, file = "output/scCATCH.csv", col.names = TRUE)