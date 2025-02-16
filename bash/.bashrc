## make a seperate R library for each version of R
export R_LIBS=`Rscript -e 'cat(file.path(Sys.getenv("HOME"), "Rlibs",paste(R.Version()$major, R.Version()$minor, sep=".")))'`
mkdir -p "$R_LIBS"

alias em='emacs -nw'
alias yes='echo "Dont be so punchy!"'
alias lj='ls -lghS'
alias duh='du -h -d 1'
alias lem='emacs -nw -q --load ~/.emacs.d/lite.el'
alias clean='rm *~ \#* 2> /dev/null'
alias dush='du -h --max-depth=1 . | sort -h -k1,1'

rem() {
    # Check if the configuration file exists  
    if [ ! -f ~/.emacs.d/lite.el ]; then  
        echo "Warning: ~/.emacs.d/lite.el not found. Emacs will run without it."  
    fi
    
    emacs "$1" -nw -q --load ~/.emacs.d/lite.el --eval '(setq buffer-read-only t)'                                                                                      
}

set_git_ssh_remote() {
    # Check if the current directory is a Git repository
    if [ ! -d ".git" ]; then
        echo "Error: This is not a Git repository."
        return 1
    fi

    # Get the folder name (assumed to be the repo name)
    repo_name=$(basename "$PWD")

    # Set the Git remote URL to SSH
    git remote set-url origin git@github.com:walaj/$repo_name.git

    # Verify the change
    git remote -v
}

