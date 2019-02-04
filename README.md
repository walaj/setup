
### Mac
Use default Terminal application. Set "option as meta" with Terminal > Preferences > Profiles > Basic (or whatever profile) > Keyboard > Turn on Option as Meta

### Emacs

Make sure Emacs is version 26 or greater
``emacs --version``

Make the ``emacs.d`` directory in home
``mkdir -p ~/.emacs.d``

Move the file ``emacs/init.d`` from this repos into ``~/.emacs.d``. The first time it opens, it should compile and install several packages. See ``emacs`` Github folder for more directions.

### Tmux
