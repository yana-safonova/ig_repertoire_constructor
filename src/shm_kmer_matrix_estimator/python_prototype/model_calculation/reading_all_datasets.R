# Loading datasets
## Loading GMC

setwd("~/chihua")

shm.k.mer.model.reader <- function(filepath) {
  model <- read.csv(filepath)
  rownames(model) <- model$X
  model$X <- NULL
  #cbind(model / rowSums(model), N=rowSums(model))
  model
}

# yale.shm <- read.csv('/home/ashlemov/Yandex.Disk/Documents/lab/vers2/src/shm_kmer_model/output/yale_substitution.csv', sep=" ")
# rownames(yale.shm) <- yale.shm$Fivemer
# yale.shm$Fivemer <- NULL
# head(yale.shm)

# syn <- which(yale.shm$Source == "Measured")

GMC.1.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/1/shm_model_utils/model_trivial_mm.csv")
GMC.2.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/2/shm_model_utils/model_trivial_mm.csv")
GMC.3.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/3/shm_model_utils/model_trivial_mm.csv")
GMC.4.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/4/shm_model_utils/model_trivial_mm.csv")
GMC.5.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/5/shm_model_utils/model_trivial_mm.csv")
GMC.6.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/6/shm_model_utils/model_trivial_mm.csv")
GMC.7.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/7/shm_model_utils/model_trivial_mm.csv")
GMC.8.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/8/shm_model_utils/model_trivial_mm.csv")
GMC.9.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/9/shm_model_utils/model_trivial_mm.csv")
GMC.10.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/10/shm_model_utils/model_trivial_mm.csv")

GMC.IG.trivial <- list(
  GMC.1.IG.trivial,
  GMC.2.IG.trivial,
  GMC.3.IG.trivial,
  GMC.4.IG.trivial,
  GMC.5.IG.trivial,
  GMC.6.IG.trivial,
  GMC.7.IG.trivial,
  GMC.8.IG.trivial,
  GMC.9.IG.trivial,
  GMC.10.IG.trivial
)



GMC.1.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/1/shm_model_utils/model_k_neighbour_mm.csv")
GMC.2.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/2/shm_model_utils/model_k_neighbour_mm.csv")
GMC.3.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/3/shm_model_utils/model_k_neighbour_mm.csv")
GMC.4.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/4/shm_model_utils/model_k_neighbour_mm.csv")
GMC.5.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/5/shm_model_utils/model_k_neighbour_mm.csv")
GMC.6.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/6/shm_model_utils/model_k_neighbour_mm.csv")
GMC.7.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/7/shm_model_utils/model_k_neighbour_mm.csv")
GMC.8.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/8/shm_model_utils/model_k_neighbour_mm.csv")
GMC.9.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/9/shm_model_utils/model_k_neighbour_mm.csv")
GMC.10.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IG/10/shm_model_utils/model_k_neighbour_mm.csv")

GMC.IG.k_neighbour <- list(
  GMC.1.IG.k_neighbour,
  GMC.2.IG.k_neighbour,
  GMC.3.IG.k_neighbour,
  GMC.4.IG.k_neighbour,
  GMC.5.IG.k_neighbour,
  GMC.6.IG.k_neighbour,
  GMC.7.IG.k_neighbour,
  GMC.8.IG.k_neighbour,
  GMC.9.IG.k_neighbour,
  GMC.10.IG.k_neighbour
)



GMC.1.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/1/shm_model_utils/model_trivial_mm.csv")
GMC.2.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/2/shm_model_utils/model_trivial_mm.csv")
GMC.3.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/3/shm_model_utils/model_trivial_mm.csv")
GMC.4.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/4/shm_model_utils/model_trivial_mm.csv")
GMC.5.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/5/shm_model_utils/model_trivial_mm.csv")
GMC.6.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/6/shm_model_utils/model_trivial_mm.csv")
GMC.7.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/7/shm_model_utils/model_trivial_mm.csv")
GMC.8.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/8/shm_model_utils/model_trivial_mm.csv")
GMC.9.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/9/shm_model_utils/model_trivial_mm.csv")
GMC.10.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/10/shm_model_utils/model_trivial_mm.csv")

GMC.IGK.trivial <- list(
  GMC.1.IGK.trivial,
  GMC.2.IGK.trivial,
  GMC.3.IGK.trivial,
  GMC.4.IGK.trivial,
  GMC.5.IGK.trivial,
  GMC.6.IGK.trivial,
  GMC.7.IGK.trivial,
  GMC.8.IGK.trivial,
  GMC.9.IGK.trivial,
  GMC.10.IGK.trivial
)



GMC.1.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/1/shm_model_utils/model_k_neighbour_mm.csv")
GMC.2.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/2/shm_model_utils/model_k_neighbour_mm.csv")
GMC.3.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/3/shm_model_utils/model_k_neighbour_mm.csv")
GMC.4.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/4/shm_model_utils/model_k_neighbour_mm.csv")
GMC.5.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/5/shm_model_utils/model_k_neighbour_mm.csv")
GMC.6.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/6/shm_model_utils/model_k_neighbour_mm.csv")
GMC.7.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/7/shm_model_utils/model_k_neighbour_mm.csv")
GMC.8.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/8/shm_model_utils/model_k_neighbour_mm.csv")
GMC.9.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/9/shm_model_utils/model_k_neighbour_mm.csv")
GMC.10.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGK/10/shm_model_utils/model_k_neighbour_mm.csv")

GMC.IGK.k_neighbour <- list(
  GMC.1.IGK.k_neighbour,
  GMC.2.IGK.k_neighbour,
  GMC.3.IGK.k_neighbour,
  GMC.4.IGK.k_neighbour,
  GMC.5.IGK.k_neighbour,
  GMC.6.IGK.k_neighbour,
  GMC.7.IGK.k_neighbour,
  GMC.8.IGK.k_neighbour,
  GMC.9.IGK.k_neighbour,
  GMC.10.IGK.k_neighbour
)



GMC.1.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/1/shm_model_utils/model_trivial_mm.csv")
GMC.2.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/2/shm_model_utils/model_trivial_mm.csv")
GMC.3.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/3/shm_model_utils/model_trivial_mm.csv")
GMC.4.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/4/shm_model_utils/model_trivial_mm.csv")
GMC.5.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/5/shm_model_utils/model_trivial_mm.csv")
GMC.6.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/6/shm_model_utils/model_trivial_mm.csv")
GMC.7.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/7/shm_model_utils/model_trivial_mm.csv")
GMC.8.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/8/shm_model_utils/model_trivial_mm.csv")
GMC.9.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/9/shm_model_utils/model_trivial_mm.csv")
GMC.10.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/10/shm_model_utils/model_trivial_mm.csv")

GMC.IGL.trivial <- list(
  GMC.1.IGL.trivial,
  GMC.2.IGL.trivial,
  GMC.3.IGL.trivial,
  GMC.4.IGL.trivial,
  GMC.5.IGL.trivial,
  GMC.6.IGL.trivial,
  GMC.7.IGL.trivial,
  GMC.8.IGL.trivial,
  GMC.9.IGL.trivial,
  GMC.10.IGL.trivial
)



GMC.1.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/1/shm_model_utils/model_k_neighbour_mm.csv")
GMC.2.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/2/shm_model_utils/model_k_neighbour_mm.csv")
GMC.3.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/3/shm_model_utils/model_k_neighbour_mm.csv")
GMC.4.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/4/shm_model_utils/model_k_neighbour_mm.csv")
GMC.5.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/5/shm_model_utils/model_k_neighbour_mm.csv")
GMC.6.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/6/shm_model_utils/model_k_neighbour_mm.csv")
GMC.7.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/7/shm_model_utils/model_k_neighbour_mm.csv")
GMC.8.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/8/shm_model_utils/model_k_neighbour_mm.csv")
GMC.9.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/9/shm_model_utils/model_k_neighbour_mm.csv")
GMC.10.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGL/10/shm_model_utils/model_k_neighbour_mm.csv")

GMC.IGL.k_neighbour <- list(
  GMC.1.IGL.k_neighbour,
  GMC.2.IGL.k_neighbour,
  GMC.3.IGL.k_neighbour,
  GMC.4.IGL.k_neighbour,
  GMC.5.IGL.k_neighbour,
  GMC.6.IGL.k_neighbour,
  GMC.7.IGL.k_neighbour,
  GMC.8.IGL.k_neighbour,
  GMC.9.IGL.k_neighbour,
  GMC.10.IGL.k_neighbour
)



GMC.1.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/1/shm_model_utils/model_trivial_mm.csv")
GMC.2.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/2/shm_model_utils/model_trivial_mm.csv")
GMC.3.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/3/shm_model_utils/model_trivial_mm.csv")
GMC.4.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/4/shm_model_utils/model_trivial_mm.csv")
GMC.5.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/5/shm_model_utils/model_trivial_mm.csv")
GMC.6.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/6/shm_model_utils/model_trivial_mm.csv")
GMC.7.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/7/shm_model_utils/model_trivial_mm.csv")
GMC.8.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/8/shm_model_utils/model_trivial_mm.csv")
GMC.9.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/9/shm_model_utils/model_trivial_mm.csv")
GMC.10.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/10/shm_model_utils/model_trivial_mm.csv")

GMC.IGH.trivial <- list(
  GMC.1.IGH.trivial,
  GMC.2.IGH.trivial,
  GMC.3.IGH.trivial,
  GMC.4.IGH.trivial,
  GMC.5.IGH.trivial,
  GMC.6.IGH.trivial,
  GMC.7.IGH.trivial,
  GMC.8.IGH.trivial,
  GMC.9.IGH.trivial,
  GMC.10.IGH.trivial
)



GMC.1.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/1/shm_model_utils/model_k_neighbour_mm.csv")
GMC.2.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/2/shm_model_utils/model_k_neighbour_mm.csv")
GMC.3.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/3/shm_model_utils/model_k_neighbour_mm.csv")
GMC.4.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/4/shm_model_utils/model_k_neighbour_mm.csv")
GMC.5.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/5/shm_model_utils/model_k_neighbour_mm.csv")
GMC.6.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/6/shm_model_utils/model_k_neighbour_mm.csv")
GMC.7.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/7/shm_model_utils/model_k_neighbour_mm.csv")
GMC.8.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/8/shm_model_utils/model_k_neighbour_mm.csv")
GMC.9.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/9/shm_model_utils/model_k_neighbour_mm.csv")
GMC.10.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/GMC/VJF_output_loci_IGH/10/shm_model_utils/model_k_neighbour_mm.csv")

GMC.IGH.k_neighbour <- list(
  GMC.1.IGH.k_neighbour,
  GMC.2.IGH.k_neighbour,
  GMC.3.IGH.k_neighbour,
  GMC.4.IGH.k_neighbour,
  GMC.5.IGH.k_neighbour,
  GMC.6.IGH.k_neighbour,
  GMC.7.IGH.k_neighbour,
  GMC.8.IGH.k_neighbour,
  GMC.9.IGH.k_neighbour,
  GMC.10.IGH.k_neighbour
)


## Loading IDO


IDO.11.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/11/shm_model_utils/model_trivial_mm.csv")
IDO.12.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/12/shm_model_utils/model_trivial_mm.csv")
IDO.13.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/13/shm_model_utils/model_trivial_mm.csv")
IDO.14.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/14/shm_model_utils/model_trivial_mm.csv")
IDO.15.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/15/shm_model_utils/model_trivial_mm.csv")
IDO.16.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/16/shm_model_utils/model_trivial_mm.csv")
IDO.17.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/17/shm_model_utils/model_trivial_mm.csv")
IDO.18.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/18/shm_model_utils/model_trivial_mm.csv")
IDO.19.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/19/shm_model_utils/model_trivial_mm.csv")
IDO.20.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/20/shm_model_utils/model_trivial_mm.csv")

IDO.IG.trivial <- list(
  IDO.11.IG.trivial,
  IDO.12.IG.trivial,
  IDO.13.IG.trivial,
  IDO.14.IG.trivial,
  IDO.15.IG.trivial,
  IDO.16.IG.trivial,
  IDO.17.IG.trivial,
  IDO.18.IG.trivial,
  IDO.19.IG.trivial,
  IDO.20.IG.trivial
)



IDO.11.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/11/shm_model_utils/model_k_neighbour_mm.csv")
IDO.12.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/12/shm_model_utils/model_k_neighbour_mm.csv")
IDO.13.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/13/shm_model_utils/model_k_neighbour_mm.csv")
IDO.14.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/14/shm_model_utils/model_k_neighbour_mm.csv")
IDO.15.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/15/shm_model_utils/model_k_neighbour_mm.csv")
IDO.16.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/16/shm_model_utils/model_k_neighbour_mm.csv")
IDO.17.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/17/shm_model_utils/model_k_neighbour_mm.csv")
IDO.18.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/18/shm_model_utils/model_k_neighbour_mm.csv")
IDO.19.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/19/shm_model_utils/model_k_neighbour_mm.csv")
IDO.20.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IG/20/shm_model_utils/model_k_neighbour_mm.csv")

IDO.IG.k_neighbour <- list(
  IDO.11.IG.k_neighbour,
  IDO.12.IG.k_neighbour,
  IDO.13.IG.k_neighbour,
  IDO.14.IG.k_neighbour,
  IDO.15.IG.k_neighbour,
  IDO.16.IG.k_neighbour,
  IDO.17.IG.k_neighbour,
  IDO.18.IG.k_neighbour,
  IDO.19.IG.k_neighbour,
  IDO.20.IG.k_neighbour
)



IDO.11.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/11/shm_model_utils/model_trivial_mm.csv")
IDO.12.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/12/shm_model_utils/model_trivial_mm.csv")
IDO.13.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/13/shm_model_utils/model_trivial_mm.csv")
IDO.14.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/14/shm_model_utils/model_trivial_mm.csv")
IDO.15.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/15/shm_model_utils/model_trivial_mm.csv")
IDO.16.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/16/shm_model_utils/model_trivial_mm.csv")
IDO.17.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/17/shm_model_utils/model_trivial_mm.csv")
IDO.18.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/18/shm_model_utils/model_trivial_mm.csv")
IDO.19.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/19/shm_model_utils/model_trivial_mm.csv")
IDO.20.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/20/shm_model_utils/model_trivial_mm.csv")

IDO.IGK.trivial <- list(
  IDO.11.IGK.trivial,
  IDO.12.IGK.trivial,
  IDO.13.IGK.trivial,
  IDO.14.IGK.trivial,
  IDO.15.IGK.trivial,
  IDO.16.IGK.trivial,
  IDO.17.IGK.trivial,
  IDO.18.IGK.trivial,
  IDO.19.IGK.trivial,
  IDO.20.IGK.trivial
)



IDO.11.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/11/shm_model_utils/model_k_neighbour_mm.csv")
IDO.12.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/12/shm_model_utils/model_k_neighbour_mm.csv")
IDO.13.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/13/shm_model_utils/model_k_neighbour_mm.csv")
IDO.14.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/14/shm_model_utils/model_k_neighbour_mm.csv")
IDO.15.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/15/shm_model_utils/model_k_neighbour_mm.csv")
IDO.16.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/16/shm_model_utils/model_k_neighbour_mm.csv")
IDO.17.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/17/shm_model_utils/model_k_neighbour_mm.csv")
IDO.18.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/18/shm_model_utils/model_k_neighbour_mm.csv")
IDO.19.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/19/shm_model_utils/model_k_neighbour_mm.csv")
IDO.20.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGK/20/shm_model_utils/model_k_neighbour_mm.csv")

IDO.IGK.k_neighbour <- list(
  IDO.11.IGK.k_neighbour,
  IDO.12.IGK.k_neighbour,
  IDO.13.IGK.k_neighbour,
  IDO.14.IGK.k_neighbour,
  IDO.15.IGK.k_neighbour,
  IDO.16.IGK.k_neighbour,
  IDO.17.IGK.k_neighbour,
  IDO.18.IGK.k_neighbour,
  IDO.19.IGK.k_neighbour,
  IDO.20.IGK.k_neighbour
)



IDO.11.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/11/shm_model_utils/model_trivial_mm.csv")
IDO.12.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/12/shm_model_utils/model_trivial_mm.csv")
IDO.13.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/13/shm_model_utils/model_trivial_mm.csv")
IDO.14.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/14/shm_model_utils/model_trivial_mm.csv")
IDO.15.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/15/shm_model_utils/model_trivial_mm.csv")
IDO.16.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/16/shm_model_utils/model_trivial_mm.csv")
IDO.17.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/17/shm_model_utils/model_trivial_mm.csv")
IDO.18.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/18/shm_model_utils/model_trivial_mm.csv")
IDO.19.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/19/shm_model_utils/model_trivial_mm.csv")
IDO.20.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/20/shm_model_utils/model_trivial_mm.csv")


IDO.IGL.trivial <- list(
  IDO.11.IGL.trivial,
  IDO.12.IGL.trivial,
  IDO.13.IGL.trivial,
  IDO.14.IGL.trivial,
  IDO.15.IGL.trivial,
  IDO.16.IGL.trivial,
  IDO.17.IGL.trivial,
  IDO.18.IGL.trivial,
  IDO.19.IGL.trivial,
  IDO.20.IGL.trivial
)



IDO.11.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/11/shm_model_utils/model_k_neighbour_mm.csv")
IDO.12.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/12/shm_model_utils/model_k_neighbour_mm.csv")
IDO.13.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/13/shm_model_utils/model_k_neighbour_mm.csv")
IDO.14.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/14/shm_model_utils/model_k_neighbour_mm.csv")
IDO.15.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/15/shm_model_utils/model_k_neighbour_mm.csv")
IDO.16.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/16/shm_model_utils/model_k_neighbour_mm.csv")
IDO.17.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/17/shm_model_utils/model_k_neighbour_mm.csv")
IDO.18.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/18/shm_model_utils/model_k_neighbour_mm.csv")
IDO.19.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/19/shm_model_utils/model_k_neighbour_mm.csv")
IDO.20.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGL/20/shm_model_utils/model_k_neighbour_mm.csv")

IDO.IGL.k_neighbour <- list(
  IDO.11.IGL.k_neighbour,
  IDO.12.IGL.k_neighbour,
  IDO.13.IGL.k_neighbour,
  IDO.14.IGL.k_neighbour,
  IDO.15.IGL.k_neighbour,
  IDO.16.IGL.k_neighbour,
  IDO.17.IGL.k_neighbour,
  IDO.18.IGL.k_neighbour,
  IDO.19.IGL.k_neighbour,
  IDO.20.IGL.k_neighbour
)



IDO.11.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/11/shm_model_utils/model_trivial_mm.csv")
IDO.12.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/12/shm_model_utils/model_trivial_mm.csv")
IDO.13.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/13/shm_model_utils/model_trivial_mm.csv")
IDO.14.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/14/shm_model_utils/model_trivial_mm.csv")
IDO.15.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/15/shm_model_utils/model_trivial_mm.csv")
IDO.16.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/16/shm_model_utils/model_trivial_mm.csv")
IDO.17.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/17/shm_model_utils/model_trivial_mm.csv")
IDO.18.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/18/shm_model_utils/model_trivial_mm.csv")
IDO.19.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/19/shm_model_utils/model_trivial_mm.csv")
IDO.20.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/20/shm_model_utils/model_trivial_mm.csv")


IDO.IGH.trivial <- list(
  IDO.11.IGH.trivial,
  IDO.12.IGH.trivial,
  IDO.13.IGH.trivial,
  IDO.14.IGH.trivial,
  IDO.15.IGH.trivial,
  IDO.16.IGH.trivial,
  IDO.17.IGH.trivial,
  IDO.18.IGH.trivial,
  IDO.19.IGH.trivial,
  IDO.20.IGH.trivial
)



IDO.11.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/11/shm_model_utils/model_k_neighbour_mm.csv")
IDO.12.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/12/shm_model_utils/model_k_neighbour_mm.csv")
IDO.13.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/13/shm_model_utils/model_k_neighbour_mm.csv")
IDO.14.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/14/shm_model_utils/model_k_neighbour_mm.csv")
IDO.15.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/15/shm_model_utils/model_k_neighbour_mm.csv")
IDO.16.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/16/shm_model_utils/model_k_neighbour_mm.csv")
IDO.17.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/17/shm_model_utils/model_k_neighbour_mm.csv")
IDO.18.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/18/shm_model_utils/model_k_neighbour_mm.csv")
IDO.19.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/19/shm_model_utils/model_k_neighbour_mm.csv")
IDO.20.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/IDO/VJF_output_loci_IGH/20/shm_model_utils/model_k_neighbour_mm.csv")

IDO.IGH.k_neighbour <- list(
  IDO.11.IGH.k_neighbour,
  IDO.12.IGH.k_neighbour,
  IDO.13.IGH.k_neighbour,
  IDO.14.IGH.k_neighbour,
  IDO.15.IGH.k_neighbour,
  IDO.16.IGH.k_neighbour,
  IDO.17.IGH.k_neighbour,
  IDO.18.IGH.k_neighbour,
  IDO.19.IGH.k_neighbour,
  IDO.20.IGH.k_neighbour
)


## Loading FV 


FV.21.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/21/shm_model_utils/model_trivial_mm.csv")
FV.22.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/22/shm_model_utils/model_trivial_mm.csv")
FV.23.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/23/shm_model_utils/model_trivial_mm.csv")
FV.24.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/24/shm_model_utils/model_trivial_mm.csv")
FV.25.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/25/shm_model_utils/model_trivial_mm.csv")
FV.26.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/26/shm_model_utils/model_trivial_mm.csv")
FV.27.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/27/shm_model_utils/model_trivial_mm.csv")
FV.28.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/28/shm_model_utils/model_trivial_mm.csv")
FV.29.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/29/shm_model_utils/model_trivial_mm.csv")
FV.30.IG.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/30/shm_model_utils/model_trivial_mm.csv")

FV.IG.trivial <- list(
  FV.21.IG.trivial,
  FV.22.IG.trivial,
  FV.23.IG.trivial,
  FV.24.IG.trivial,
  FV.25.IG.trivial,
  FV.26.IG.trivial,
  FV.27.IG.trivial,
  FV.28.IG.trivial,
  FV.29.IG.trivial,
  FV.30.IG.trivial
)



FV.21.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/21/shm_model_utils/model_k_neighbour_mm.csv")
FV.22.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/22/shm_model_utils/model_k_neighbour_mm.csv")
FV.23.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/23/shm_model_utils/model_k_neighbour_mm.csv")
FV.24.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/24/shm_model_utils/model_k_neighbour_mm.csv")
FV.25.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/25/shm_model_utils/model_k_neighbour_mm.csv")
FV.26.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/26/shm_model_utils/model_k_neighbour_mm.csv")
FV.27.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/27/shm_model_utils/model_k_neighbour_mm.csv")
FV.28.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/28/shm_model_utils/model_k_neighbour_mm.csv")
FV.29.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/29/shm_model_utils/model_k_neighbour_mm.csv")
FV.30.IG.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IG/30/shm_model_utils/model_k_neighbour_mm.csv")

FV.IG.k_neighbour <- list(
  FV.21.IG.k_neighbour,
  FV.22.IG.k_neighbour,
  FV.23.IG.k_neighbour,
  FV.24.IG.k_neighbour,
  FV.25.IG.k_neighbour,
  FV.26.IG.k_neighbour,
  FV.27.IG.k_neighbour,
  FV.28.IG.k_neighbour,
  FV.29.IG.k_neighbour,
  FV.30.IG.k_neighbour
)



FV.21.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/21/shm_model_utils/model_trivial_mm.csv")
FV.22.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/22/shm_model_utils/model_trivial_mm.csv")
FV.23.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/23/shm_model_utils/model_trivial_mm.csv")
FV.24.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/24/shm_model_utils/model_trivial_mm.csv")
FV.25.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/25/shm_model_utils/model_trivial_mm.csv")
FV.26.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/26/shm_model_utils/model_trivial_mm.csv")
FV.27.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/27/shm_model_utils/model_trivial_mm.csv")
FV.28.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/28/shm_model_utils/model_trivial_mm.csv")
FV.29.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/29/shm_model_utils/model_trivial_mm.csv")
FV.30.IGK.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/30/shm_model_utils/model_trivial_mm.csv")

FV.IGK.trivial <- list(
  FV.21.IGK.trivial,
  FV.22.IGK.trivial,
  FV.23.IGK.trivial,
  FV.24.IGK.trivial,
  FV.25.IGK.trivial,
  FV.26.IGK.trivial,
  FV.27.IGK.trivial,
  FV.28.IGK.trivial,
  FV.29.IGK.trivial,
  FV.30.IGK.trivial
)



FV.21.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/21/shm_model_utils/model_k_neighbour_mm.csv")
FV.22.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/22/shm_model_utils/model_k_neighbour_mm.csv")
FV.23.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/23/shm_model_utils/model_k_neighbour_mm.csv")
FV.24.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/24/shm_model_utils/model_k_neighbour_mm.csv")
FV.25.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/25/shm_model_utils/model_k_neighbour_mm.csv")
FV.26.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/26/shm_model_utils/model_k_neighbour_mm.csv")
FV.27.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/27/shm_model_utils/model_k_neighbour_mm.csv")
FV.28.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/28/shm_model_utils/model_k_neighbour_mm.csv")
FV.29.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/29/shm_model_utils/model_k_neighbour_mm.csv")
FV.30.IGK.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGK/30/shm_model_utils/model_k_neighbour_mm.csv")

FV.IGK.k_neighbour <- list(
  FV.21.IGK.k_neighbour,
  FV.22.IGK.k_neighbour,
  FV.23.IGK.k_neighbour,
  FV.24.IGK.k_neighbour,
  FV.25.IGK.k_neighbour,
  FV.26.IGK.k_neighbour,
  FV.27.IGK.k_neighbour,
  FV.28.IGK.k_neighbour,
  FV.29.IGK.k_neighbour,
  FV.30.IGK.k_neighbour
)



FV.21.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/21/shm_model_utils/model_trivial_mm.csv")
FV.22.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/22/shm_model_utils/model_trivial_mm.csv")
FV.23.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/23/shm_model_utils/model_trivial_mm.csv")
FV.24.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/24/shm_model_utils/model_trivial_mm.csv")
FV.25.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/25/shm_model_utils/model_trivial_mm.csv")
FV.26.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/26/shm_model_utils/model_trivial_mm.csv")
FV.27.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/27/shm_model_utils/model_trivial_mm.csv")
FV.28.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/28/shm_model_utils/model_trivial_mm.csv")
FV.29.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/29/shm_model_utils/model_trivial_mm.csv")
FV.30.IGL.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/30/shm_model_utils/model_trivial_mm.csv")

FV.IGL.trivial <- list(
  FV.21.IGL.trivial,
  FV.22.IGL.trivial,
  FV.23.IGL.trivial,
  FV.24.IGL.trivial,
  FV.25.IGL.trivial,
  FV.26.IGL.trivial,
  FV.27.IGL.trivial,
  FV.28.IGL.trivial,
  FV.29.IGL.trivial,
  FV.30.IGL.trivial
)



FV.21.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/21/shm_model_utils/model_k_neighbour_mm.csv")
FV.22.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/22/shm_model_utils/model_k_neighbour_mm.csv")
FV.23.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/23/shm_model_utils/model_k_neighbour_mm.csv")
FV.24.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/24/shm_model_utils/model_k_neighbour_mm.csv")
FV.25.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/25/shm_model_utils/model_k_neighbour_mm.csv")
FV.26.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/26/shm_model_utils/model_k_neighbour_mm.csv")
FV.27.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/27/shm_model_utils/model_k_neighbour_mm.csv")
FV.28.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/28/shm_model_utils/model_k_neighbour_mm.csv")
FV.29.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/29/shm_model_utils/model_k_neighbour_mm.csv")
FV.30.IGL.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGL/30/shm_model_utils/model_k_neighbour_mm.csv")

FV.IGL.k_neighbour <- list(
  FV.21.IGL.k_neighbour,
  FV.22.IGL.k_neighbour,
  FV.23.IGL.k_neighbour,
  FV.24.IGL.k_neighbour,
  FV.25.IGL.k_neighbour,
  FV.26.IGL.k_neighbour,
  FV.27.IGL.k_neighbour,
  FV.28.IGL.k_neighbour,
  FV.29.IGL.k_neighbour,
  FV.30.IGL.k_neighbour
)



FV.21.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/21/shm_model_utils/model_trivial_mm.csv")
FV.22.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/22/shm_model_utils/model_trivial_mm.csv")
FV.23.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/23/shm_model_utils/model_trivial_mm.csv")
FV.24.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/24/shm_model_utils/model_trivial_mm.csv")
FV.25.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/25/shm_model_utils/model_trivial_mm.csv")
FV.26.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/26/shm_model_utils/model_trivial_mm.csv")
FV.27.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/27/shm_model_utils/model_trivial_mm.csv")
FV.28.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/28/shm_model_utils/model_trivial_mm.csv")
FV.29.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/29/shm_model_utils/model_trivial_mm.csv")
FV.30.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/30/shm_model_utils/model_trivial_mm.csv")

FV.IGH.trivial <- list(
  FV.21.IGH.trivial,
  FV.22.IGH.trivial,
  FV.23.IGH.trivial,
  FV.24.IGH.trivial,
  FV.25.IGH.trivial,
  FV.26.IGH.trivial,
  FV.27.IGH.trivial,
  FV.28.IGH.trivial,
  FV.29.IGH.trivial,
  FV.30.IGH.trivial
)



FV.21.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/21/shm_model_utils/model_k_neighbour_mm.csv")
FV.22.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/22/shm_model_utils/model_k_neighbour_mm.csv")
FV.23.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/23/shm_model_utils/model_k_neighbour_mm.csv")
FV.24.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/24/shm_model_utils/model_k_neighbour_mm.csv")
FV.25.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/25/shm_model_utils/model_k_neighbour_mm.csv")
FV.26.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/26/shm_model_utils/model_k_neighbour_mm.csv")
FV.27.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/27/shm_model_utils/model_k_neighbour_mm.csv")
FV.28.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/28/shm_model_utils/model_k_neighbour_mm.csv")
FV.29.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/29/shm_model_utils/model_k_neighbour_mm.csv")
FV.30.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/AbVitro/flu_time_course/FV/VJF_output_loci_IGH/30/shm_model_utils/model_k_neighbour_mm.csv")

FV.IGH.k_neighbour <- list(
  FV.21.IGH.k_neighbour,
  FV.22.IGH.k_neighbour,
  FV.23.IGH.k_neighbour,
  FV.24.IGH.k_neighbour,
  FV.25.IGH.k_neighbour,
  FV.26.IGH.k_neighbour,
  FV.27.IGH.k_neighbour,
  FV.28.IGH.k_neighbour,
  FV.29.IGH.k_neighbour,
  FV.30.IGH.k_neighbour
)


age.1.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/1/shm_model_utils/model_k_neighbour_mm.csv")
age.2.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/2/shm_model_utils/model_k_neighbour_mm.csv")
age.3.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/3/shm_model_utils/model_k_neighbour_mm.csv")
age.4.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/4/shm_model_utils/model_k_neighbour_mm.csv")
age.5.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/5/shm_model_utils/model_k_neighbour_mm.csv")
age.6.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/6/shm_model_utils/model_k_neighbour_mm.csv")
age.7.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/7/shm_model_utils/model_k_neighbour_mm.csv")
age.8.IGH.k_neighbour <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/8/shm_model_utils/model_k_neighbour_mm.csv")

age.IGH.k_neighbour <- list(
  age.1.IGH.k_neighbour,
  age.2.IGH.k_neighbour,
  age.3.IGH.k_neighbour,
  age.4.IGH.k_neighbour,
  age.5.IGH.k_neighbour,
  age.6.IGH.k_neighbour,
  age.7.IGH.k_neighbour,
  age.8.IGH.k_neighbour
)

age.1.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/1/shm_model_utils/model_trivial_mm.csv")
age.2.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/2/shm_model_utils/model_trivial_mm.csv")
age.3.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/3/shm_model_utils/model_trivial_mm.csv")
age.4.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/4/shm_model_utils/model_trivial_mm.csv")
age.5.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/5/shm_model_utils/model_trivial_mm.csv")
age.6.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/6/shm_model_utils/model_trivial_mm.csv")
age.7.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/7/shm_model_utils/model_trivial_mm.csv")
age.8.IGH.trivial <- shm.k.mer.model.reader("./Sid/abzikadze/datasets/age/vjf_output_on_final_repertoire/IGH/8/shm_model_utils/model_trivial_mm.csv")

age.IGH.trivial <- list(
  age.1.IGH.trivial,
  age.2.IGH.trivial,
  age.3.IGH.trivial,
  age.4.IGH.trivial,
  age.5.IGH.trivial,
  age.6.IGH.trivial,
  age.7.IGH.trivial,
  age.8.IGH.trivial
)