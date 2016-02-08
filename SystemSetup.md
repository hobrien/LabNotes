# Flor
- Set up python
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
    - install [MySQLdb](http://mysql-python.sourceforge.net/MySQLdb.html): ```pip install MySQL-python```
        - I'm having problems with you scripts because queries need a tuple of parameters, meaning single parameters need to be followed by a comma
    - install [ETE2](http://etetoolkit.org): ```pip install ete2```

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
    - this seems to be working well. see Workflows/BegoniaCandidates.py for an example

- Install [Entrez Direct](http://www.ncbi.nlm.nih.gov/books/NBK179288):
    ```
    cd ~
  perl -MNet::FTP -e \
    '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1); $ftp->login;
     $ftp->binary; $ftp->get("/entrez/entrezdirect/edirect.zip");'
  unzip -u -q edirect.zip
  rm edirect.zip
  export PATH=$PATH:$HOME/edirect
  ./edirect/setup.sh
  ```
 
- Installed [Newick Utilities](http://cegg.unige.ch/newick_utils)
    - downloaded Darwin binary and copied utils to /usr/local/bin:
        - ```sudo cp src/nw_* /usr/local/bin```
    - There's something about a LibXML enabled version, but I haven't tried messing with it    
        
- Installed [TreeLink](http://www.treelinkapp.com) (downloaded and copied to applications)
    - treelink.html need to be [modified](https://github.com/allendecid/TreeLink/issues/1) to allow copy and paste of search terms
    - It is very nice to be able to work with CSVs, but it doesn't appear to be possible to change the tree root
        - nw_reroot from NewickUtils can do this, but it's still a pain
        - There is something called TreeLite that offers this functionality, as well as search and smooth zooming and scrolling, but without the CSV feature. Hopefully the rooting will migrate to TreeLink
        
- Merge Papers Libraries
    - At some point, BTsync stopped syncing between my work computer and Flor
    - I copied the Papers library from my work backup to the default location on Flor (/Users/heo3/Library/Application Support), but I need to check if there are any files in BTsync that aren't in the backed up version:
        - ```diff -rd ~/Library/Application\ Support/Papers2/ ~/BTSync/Papers2/ |grep -v -e 'DS_Store' -e 'Thumbs' -e '/: \.' | sort > diffs.txt```
        - I used grep to find only files that were only in the BTSync and all of them were due to name changes, not missing files. 
        - I guess it should be safe to delete the folder in BTSync now
        - I should check to make sure the rest of the folders in there are also up to date

- Merge other folders
    - /Volumes/Geinitz_backup/Users/HeathOBrien/Documents/BTSync -> /Users/heo3/Documents/BTSync
    - /Volumes/Geinitz_backup/Users/HeathOBrien/BTSync2 -> /Users/heo3/Documents/BTSync
    - /Volumes/Geinitz_backup/Users/HeathOBrien/BTSync2/WorkDesktop -> /Users/heo3/Desktop/WorkDesktop
    -  /Volumes/Geinitz_backup/Users/HeathOBrien/Desktop -> /Users/heo3/Desktop/WorkDesktop
    - /Volumes/Geinitz_backup/Users/HeathOBrien/BitTorrent Sync -> /Users/heo3/Documents
    - /Volumes/Geinitz_backup/Users/HeathOBrien/CLC_Data -> /Users/heo3/CLC_Data
    - Unmerged
        - Bioinformatics
        - Databases
        
