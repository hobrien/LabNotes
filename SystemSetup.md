# Flor
- Sort out the problem of scripts written on my iMac not working on my laptop
    - created a symlink so that anaconda is available in my home directory:
        - ```ln -s /anaconda ~/anaconda```
    - I also replaced /anaconda with ~/anaconda in my PATH in ~/.profile 
    - this still didn't solve the problem because the user name is different and I can't use relative paths
    - the solution turns out to be replacing the path with "#!/usr/bin/env python"
    - the script runs now, but it can't find biopython
    - it works if I use /usr/local/bin/python but not with ~/anaconda/python
    - install biopython with conda:
        - ```conda install biopython```

- Install a tool to preview GitHub markdown (the viewers I have don't render everything the same):
    - [grip](https://github.com/joeyespo/grip) appears to be pretty nice
    - launches preview in the default browser and updates automatically on save
    - installation:
        - ```pip install grip```
    - to use:
        - ```grip -b FILENAME```
        
- Figure out an easier way to build pipelines than Make
    - Make is a nice tool to build bioinformatics pipelines because it handles dependencies meaning that when changes are made to analyses or data, only downstream analyses are repeated
    - however, it's designed to handle software dependencies and is pretty cumbersome to use for pipelines
    - there is a [long list](https://github.com/pditommaso/awesome-pipeline) of alternative tools available
    - of all these, [flex](https://github.com/druths/flex) seems promising as it's billed as Make for data science
    - installation is easy:
        - ```pip install flexds```
    - this has some very cool features (primarily combining python and bash commands in an easy way), but it doesn't seem to support what I really want, which is a way to process a folder full of files without having to list each of them
    - after spending quite some time browsing that list of tools, I've settled on [doit](http://pydoit.org/)
        - this is a python package that seems to offer some of the functionality I'm looking for, but crucially, there's a [Software Carpentry](http://swcarpentry.github.io/bc/intermediate/doit) tutorial for it, which is the way I learned Make, so I'm hoping it will cover similar ground
        - install:
            ```pip install doit```
                    