;; turn off the menu bar
(menu-bar-mode -1)

(set-keyboard-coding-system nil)

;; set the toggle-truncate-lines
(set-default 'truncate-lines t)

;; set python
(setq py-python-command "/usr/local/bin/python")

;; turn on column numbering
(setq column-number-mode t)

;; set a comment + time-stamp function
(defun insert-current-date ()                                                                                                                                                                                                                                                                                               
  (interactive)                                                                                                                                                                                                                                                                                                            
  (insert (format-time-string "#' %A, %b %d, %Y %r")))                                                                                                                                                                                                                                                                     
(global-set-key (kbd "C-c d") `insert-current-date)

(setq scroll-step 1 scroll-conservatively 10000)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; vi like parenthesis jump
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
;(global-set-key "M-9" 'kill-whole-line)   ; similar to vim dd
(define-key ctl-x-map "t" 'transpose-lines)     ; Ctrl-x t runs the
                                                ; transpose-lines function.
 
(setq require-final-newline t)                  ; Make sure file always
                                                ; ends with a newline.
;(x-set-font "-adobe-courier-medium-r-*-*-20-*-*-*-*-*-*-*")

;(setq default-major-mode 'text-mode)            ; Default mode is text-mode.

(setq text-mode-hook                            ; Enable auto-fill-mode
 '(lambda () (auto-fill-mode 1)))               ; whenever using text-mode.

(setq delete-auto-save-files t)                 ; Delete unnecessary


; Marcin stuff
(set-background-color "grey20")
(set-foreground-color "white")
(set-cursor-color "white")
(set-frame-font "-*-fixed-*-*-*--11-*-*-*-c-*-*-*")
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

(require 'package)

(add-to-list 'package-archives
       '("melpa" . "http://melpa.org/packages/") t)

(package-initialize)
(when (not package-archive-contents)
  (package-refresh-contents))

(defvar myPackages
  '(better-defaults
    elpy
    material-theme))

(mapc #'(lambda (package)
    (unless (package-installed-p package)
      (package-install package)))
      myPackages)
(custom-set-variables
 ;; custom-set-variables was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 '(package-selected-packages (quote (material-theme better-defaults))))
(custom-set-faces
 ;; custom-set-faces was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 )

;; load the materal theme
;; (load-theme 'material t) ;; load material theme

;; load the pthony environment
(elpy-enable)
