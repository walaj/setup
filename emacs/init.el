(require 'package)

(setq package-archives
      '(("melpa" . "https://melpa.org/packages/")
        ("gnu"   . "https://elpa.gnu.org/packages/")))

(unless package-archive-contents
  (package-refresh-contents))

(unless (package-installed-p 'use-package)
  (package-install 'use-package))

(require 'use-package)
(setq use-package-always-ensure t)

;; turn off the menu bar
(menu-bar-mode -1)

;; go right to scratch buffer
(setq inhibit-startup-message t)

;; display line numbers
(global-display-line-numbers-mode 1)

;; Make ESC quit prompts
(global-set-key (kbd "<escape>") 'keyboard-escape-quit)

;; but disable line numbers for some modes
(dolist (mode '(org-mode-hook
		term-mode-hook
		eshell-mode-hook
		treemacs-mode-hook
		shell-mode-hook))
  (add-hook mode (lambda () (display-line-numbers-mode 0))))

;; save the minibuffer history (can scroll up like in command  line to get history)
(setq history-length 25)
(savehist-mode 1)
(save-place-mode 1)

;; place automated package customizations in file other than init.el
(setq custom-file (locate-user-emacs-file "custom-vars.el"))
(load custom-file 'noerror 'nomessage)

;; if < emacs 28 need to add in minibuffer
;M-x package-install RET modus-themes RET
;(load-theme 'modus-vivendi 1)

;; allow recent-file mode. M-x recentf-open-files
(recentf-mode 1)

;; automatically reload buffer if file changed below you 
(global-auto-revert-mode 1)

(set-keyboard-coding-system nil)

;; set the toggle-truncate-lines
(set-default 'truncate-lines t)

;; set python
;(setq py-python-command "/usr/local/bin/python")

;; turn on column numbering
(setq column-number-mode t)

;; set a comment + time-stamp function
(defun insert-current-date ()                                                                                                                                                                                                                                                                                               
  (interactive)                                                                                                                                                                                                                                                                                                            
  (insert (format-time-string "#' %A, %b %d, %Y %r")))                                                                                                                                                                                                                                                                     
(global-set-key (kbd "C-c d") `insert-current-date)

(setq scroll-step 1 scroll-conservatively 10000)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; Vi like parenthesis jump
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defun goto-match-paren (arg)
  "Go to the matching parenthesis if on parenthesis, otherwise insert %.
vi style of % jumping to matching brace."
  (interactive "p")
  (cond ((looking-at "\\s\(") (forward-list 1) (backward-char 1))
        ((looking-at "\\s\)") (forward-char 1) (backward-list 1))
        (t (self-insert-command (or arg 1)))))
(global-set-key (kbd "C-%") 'goto-match-paren)


;;;;;;; more intuitive window navigation
(global-set-key (kbd "C-x <up>") 'windmove-up)
(global-set-key (kbd "C-x <down>") 'windmove-down)
(global-set-key (kbd "C-x <right>") 'windmove-right)
(global-set-key (kbd "C-x <left>") 'windmove-left)

(transient-mark-mode t)

(define-key esc-map "G" 'goto-line)             ; Esc-G runs the goto-line
                                                ; function.
;;(global-set-key "M-9" 'kill-whole-line)   ; similar to vim dd
(define-key ctl-x-map "t" 'transpose-lines)     ; Ctrl-x t runs the
                                                ; transpose-lines function.
 
(setq require-final-newline t)                  ; Make sure file always
                                                ; ends with a newline.
;;(x-set-font "-adobe-courier-medium-r-*-*-20-*-*-*-*-*-*-*")

;;(setq default-major-mode 'text-mode)            ; Default mode is text-mode.

(setq text-mode-hook                            ; Enable auto-fill-mode
 '(lambda () (auto-fill-mode 1)))               ; whenever using text-mode.

(setq delete-auto-save-files t)                 ; Delete unnecessary


; Marcin stuff
(set-background-color "grey20")
(set-foreground-color "white")
(set-cursor-color "white")
;; (set-frame-font "-*-fixed-*-*-*--11-*-*-*-c-*-*-*")
(add-to-list 'default-frame-alist '(font . "6x13"))
(add-to-list 'default-frame-alist '(foreground-color . "white"))
(add-to-list 'default-frame-alist '(background-color . "grey20"))
(add-to-list 'default-frame-alist '(cursor-color . "white"))
(set-face-background 'default "grey20")
(set-face-foreground 'default "white")

;; ask before exiting
(defun my-exit-from-emacs()
  (interactive)
  ( if (yes-or-no-p "Do you want to exit ")
      (save-buffers-kill-emacs)))
(global-set-key "\C-x\C-c" 'my-exit-from-emacs)

;; INSTALL PACKAGES
;; --------------------------------------
(setq gnutls-algorithm-priority "NORMAL:-VERS-TLS1.3")

;; initialize package sources
(require 'package)

(setq package-archives '(("melpa" . "https://melpa.org/packages/")
			 ("org" . "https://orgmode.org/elpa/")
			 ("elpa" . "https://elpa.gnu.org/packages/")))

(package-initialize)
(unless package-archive-contents
  (package-refresh-contents))

;; install use-package if not found (M-x package-install RET ***)
;; use-package

;; use use-package to install elpy
(require 'use-package)
(setq use-package-always-ensure t) ;; install if not found

;; for integration of R into emacs
(use-package ess)

; set to desired indentation level, e.g., 4 spaces
(setq ess-indent-offset 4)
(setq ess-arg-function-offset 4)

;; Adds M-x recent command sorting for counsel-M-x
(use-package smex
:defer 1
:after counsel)

;; more rich functions in emacs
(use-package counsel
   :bind (("M-x" . counsel-M-x)
	 ("C-x b" . counsel-ibuffer)
	 ("C-x C-f" . counsel-find-file)
	 ("C-x C-b" . counsel-ibuffer)
	 :map minibuffer-local-map
	 ("C-r" . 'counsel-minibuffer-history))
 :config
 (setq ivy-initial-inputs-alist nil)) ;; Don't start searches with ^

;; more rich functions in emacs
(use-package ivy
 :bind (("C-s" . swiper)
	 :map ivy-minibuffer-map
	 ;("TAB" . ivy-alt-done)     
	 ("C-j" . ivy-next-line)
	 ("C-k" . ivy-previous-line)
	 :map ivy-switch-buffer-map
	 ("C-k" . ivy-previous-line)
	 ("C-l" . ivy-done)
	 ("C-d" . ivy-switch-buffer-kill)
	 )
 :config
 (ivy-mode 1))

;(setq user-emacs-directory "~/.cache/emacs")
;(use-package no-litering)

(use-package ivy-rich
 :init
 (ivy-rich-mode 1))

(use-package rainbow-delimiters
  :hook (prog-mode . rainbow-delimiters-mode)) 

;; display key options if delay
 (use-package which-key
  :init (which-key-mode)
  :diminish which-key-mode
  :config
  (setq which-key-idle-delay 0.5)
  (setq which-key-separator " - " ))

;;
;;(use-package doom-modeline
;;  :init (doom-modeline-mode 1))

;(use-package python-mode
;  :ensure t
;  :custom
;  (python-sxhell-interpreter "python3"))

(use-package elpy
  :ensure t
  :init
  (elpy-enable))

(add-hook 'elpy-mode-hook
    (lambda ()
    (local-unset-key (kbd "C-c C-n"))
    (define-key elpy-mode-map (kbd "C-c C-n") 'elpy-shell-send-statement-and-step)))

;; use ipython
(setq python-shell-interpreter "ipython"
      python-shell-interpreter-args "-i --simple-prompt")

;; send highlighted region to mac clipboard
(defun clip-copy ()
  (interactive)
  (when (region-active-p)
    (shell-command-on-region (region-beginning) (region-end) "pbcopy")
    (deactivate-mark)))

(add-hook 'python-mode-hook 'whitespace-mode)
(setq whitespace-style '(face tabs trailing tab-mark))

;(use-package pyvenv
;  :ensure t
;  :config
;  (pyvenv-mode t)

  ;; Set correct Python interpreter
  ;(setq pyvenv-post-activate-hooks
  ;      (list (lambda ()
  ;             (setq python-shell-interpreter (concat pyvenv-virtual-env "bin/python3")))))
  ;(setq pyvenv-post-deactivate-hooks
  ;      (list (lambda ()
  ;              (setq python-shell-interpreter "python3")))))


;; (use-package elpy;
;; OA	     :ensure t
;; 	     :init
;; 	     (elpy-enable))
