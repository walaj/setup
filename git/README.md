```
# Setup user name file
git config --global user.email "jeremiah.wala@gmail.com"
git config --global user.name "walaj"
```

Generate an RSA key (if already generated, hit "n")
``ssh-keygen -t rsa``

If it is a GitHub repository and you have administrative privileges, go to settings and click 'add SSH key'. Copy the contents of your ~/.ssh/id_rsa.pub into the field labeled 'Key'.

Don't do this with ``git config credential.helper store`` since it stores your password unencrypted

```
# Make sure repo URL is in form git+ssh, you can do this with below
git remote show origin          # check if the url is in form below
git remote set-url origin git+ssh://git@github.com/username/reponame.git
```
