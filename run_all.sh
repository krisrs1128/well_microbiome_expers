mkdir chapter/figure/pca
mkdir chapter/figure/pmd
mkdir chapter/figure/spls
mkdir chapter/figure/graph_lasso
mkdir chapter/figure/lda_cca

cd dimension_red
Rscript pca.R
Rscript pmd.R
Rscript cca.R
Rscript ccpna.R
Rscript coia.R
Rscript lda_cca.R
Rscript pca_iv.R
Rscript illustration.R
Rscript illustration_pmd.R

cd ../supervised/
Rscript spls.R
Rscript graph_lasso.R
Rscript multitask_lasso.R
