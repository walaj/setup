options(tab.width = 2)
#options(texi2dvi='/xchip/gistic/Jeremiah/software/texi2dvi/texi2dvi')
#options(texi2dvi='/usr/bin/texi2dvi4a2ps')

# set the opneing and closing functions
#.First <- function() cat("\n   Welcome to R!\n\n")
#.Last <- function()  cat("\n   Goodbye!\n\n")

## set the repos
r = getOption("repos") # hard code the UK repo for CRAN
r["CRAN"] = "http://cran.uk.r-project.org"
options(repos = r)
rm(r)

## convient function
#ucount <- function(x) {
#  length(unique(x[!is.na(x)]))
#}

#############
### functions for evaluating the total amount of memory, and per-object memory
############
# https://stackoverflow.com/questions/12391950/select-assign-to-data-table-variables-which-names-are-stored-in-a-character-ve
EVAL = function(...)eval(parse(text=paste0(...)),envir=parent.frame(2))

## per-object memory
.sizeit  <- function() format(round(sort( sapply(ls(envir=.GlobalEnv),function(x){object.size(get(x))})) / 1e6, 1), nsmall=1)

## total memory
.allsize <- function() print(object.size(x=lapply(ls(envir=.GlobalEnv), get)), units="Mb") 
