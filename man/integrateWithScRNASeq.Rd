% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integrateATACWithRNA.R
\name{integrateWithScRNASeq}
\alias{integrateWithScRNASeq}
\title{Integration with single-cell RNA sequencing data}
\usage{
integrateWithScRNASeq(
  snapSeurat,
  rnaSeurat,
  eigDims = NULL,
  snapAssay = "GeneScore",
  preprocessSnap = TRUE,
  preprocessRNA = TRUE,
  rnaTypeColnm = "ClusterName",
  reso = 0.5,
  reduction = "cca",
  nVariableFeature = 3000
)
}
\arguments{
\item{snapSeurat}{Seurat object transformed from snap object
Assume meta.data has no column named ClusterName}

\item{rnaSeurat}{Seurat obejct, reference scRNA sequencing ddta
Assume meta.data has a column for rnaType}

\item{eigDims}{dims used for snapSeurat pca,
if NULL (default), use all the dims.}

\item{snapAssay}{characters, snapSeurat assay name, default "GeneScore"}

\item{preprocessSnap}{bool, if need to normalize and scale snap, default TRUE}

\item{preprocessRNA}{bool, if need to normalize, scale rnaSeurat, default TRUE}

\item{rnaTypeColnm}{characters, colname for rnaSeurat's type, default "ClusterName"}

\item{reso}{numeric, clustering resolution parameter, default 0.5}

\item{reduction}{characters, method used in FindTransferAnchors,}

\item{nVariableFeature}{integer, default is 3000
default "cca".}
}
\value{
Seurat object, merge snapSeurat and rnaSeurat
Add predictId  and predictMaxScore from rnaSeurat for snapSeurat on metaData,
Add imputed gene expression from rnaSeurat for snapSeurat with Assay: "RNA"
Add tech as either ATAC or RNA in the meta.data
Clustering on the merged Seurat object as Co-Embedding Seurat object
}
\description{
Integration with single-cell RNA sequencing data
}
