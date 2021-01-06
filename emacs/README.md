Make a directory for Emacs settings and place ``init.el`` in it
``mkdir -p ~/.emacs.d``

See comments within ``init.d``. In particular, this sets up two systems the first time that Emacs is run after ``init.el`` is in place.

#### Python
The init.d installs ``elpy``, an Emacs package to use Python as an IDE. To use, open Python file. Run the entire file in an in-emacs REPL with ``Ctrl-c Ctrl-c``. Can also make a selction (highlight and selection with ``Ctrl-space``), then run that section with ``Ctrl-c Ctrl-c``.

#### R
Install Emacs Speak Statistics, for running R in Emacs. Start with ``M-x R``. See online tutorial for ESS for keys, but in general you can run a line with ``Ctrl-c Ctrl-n`` and a selection with ``Ctrl-c Ctrl-r``.
