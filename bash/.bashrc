## make a seperate R library for each version of R
export R_LIBS=`Rscript -e 'cat(file.path(Sys.getenv("HOME"), "Rlibs",paste(R.Version()$major, R.Version()$minor, sep=".")))'`
mkdir -p "$R_LIBS"

alias em='emacs -nw'
alias yes='echo "Dont be so punchy!"'
alias lj='ls -lghS'
