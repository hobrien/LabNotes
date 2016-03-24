# run Xcode and agree to licence terms

# this requires user input: /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew tap caskroom/cask
brew cask install firefox
brew cask install istat-menus # requires licence key eventually
brew cask install bbedit # requires licence key eventually
brew cask install lastpass
brew cask install google-chrome
brew cask install papers # requires licence key eventually
brew cask install bittorrent-sync
brew cask install github-desktop
brew cask install rstudio
brew cask install r
brew cask install vienna
brew cask install dropbox
brew cask install yojimbo # requires licence key eventually
brew cask install sequel-pro
brew cask install mactex  # need to update TexShop before it will work
brew cask install java
brew cask install skype

brew cask install anaconda
conda create -n python2 python=2.7 anaconda
# Set up bash profile
echo "~/.bash_profile" > ~/.bashrc # this doesn't appear to be working because a command using xargs couldn't find an executable in ~/bin
echo 'export PATH="~/bin:~/BTSync/Code/Python:/usr/local/sbin:/Library/TeX/texbin:$PATH"' > ~/.bash_profile
echo 'source ~/anaconda3/bin/activate python2' >> ~/.bash_profile
echo 'export PERL5LIB=/usr/local/lib/perl5/site_perl/' >> ~/.bash_profile

# fix homebrew:
sudo chown -R $(whoami) /usr/local/include
sudo chown -R $(whoami) /usr/local/lib
sudo chown -R $(whoami) /usr/local/lib/pkgconfig
sudo chown -R $(whoami) /usr/local/share/man/man5
sudo chown -R $(whoami) /usr/local/share/man/mann


brew tap homebrew/science
brew install samtools
brew install bcftools
brew install mercurial
brew install mysql
brew install pandoc
brew install igv

# Link pandoc template files
ln -s ~/BTSync/LaTEX/Templates/letter.latex /usr/local/Cellar/pandoc/1.16.0.2/share/x86_64-osx-ghc-7.10.3/pandoc-1.16.0.2/data/templates/letter.latex
ln -s ~/BTSync/LaTEX/Templates/letter.home.latex /usr/local/Cellar/pandoc/1.16.0.2/share/x86_64-osx-ghc-7.10.3/pandoc-1.16.0.2/data/templates/letter.home.latex
ln -s ~/BTSync/LaTEX/Templates/report.latex /usr/local/Cellar/pandoc/1.16.0.2/share/x86_64-osx-ghc-7.10.3/pandoc-1.16.0.2/data/templates/report.latex
# usage: pandoc -s -S -t latex --template=letter -o OUT.pdf IN.md

# install Mail ActOn
# install [checkVCF.py](https://github.com/zhanxw/checkVCF)
#    - symlink the script to ~/bin
#    - requires ~/Documents/src/

pip install grip
conda install mysql-connector-python

# this doesn't work no matter how I try to install it: pip install bx-python

