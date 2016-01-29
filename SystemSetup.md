# Flor
- Sort out the problem of scripts written on my iMac not working on my laptop
    - created a symlink so that anaconda is available in my home directory
    - ```ln -s /anaconda ~/anaconda```
    - I also replaced /anaconda with ~/anaconda in my PATH in ~/.profile 
    - this still didn't solve the problem because the user name is different and I can't use relative paths
    - the solution turns out to be replacing the path with "#!/usr/bin/env python"
    - the script runs now, but it can't find biopython
    - it works if I use /usr/local/bin/python but not with ~/anaconda/python
    