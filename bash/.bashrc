## make a seperate R library for each version of R
export R_LIBS=`Rscript -e 'cat(file.path("~/Rlibs",paste(R.Version()$major, R.Version()$minor, sep=".")))'`
