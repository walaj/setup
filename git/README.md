```
# Setup user name file
git config --global user.email "jeremiah.wala@gmail.com"
git config --global user.name "walaj"
```

Generate an RSA key (if already generated, hit "n")
``ssh-keygen -t rsa``

If it is a GitHub repository and you have administrative privileges, go to settings and click 'add SSH key'. Copy the contents of your ~/.ssh/id_rsa.pub into the field labeled 'Key'.
