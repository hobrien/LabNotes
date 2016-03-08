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

brew cask install anaconda
conda create -n python2 python=2.7 anaconda
# Set up bash profile
echo "~/.bash_profile" > ~/.bashrc # this doesn't appear to be working because a command using xargs couldn't find an executable in ~/bin
echo 'export PATH="~/bin:~/BTSync/Code/Python:/usr/local/sbin:$PATH"' > ~/.bash_profile
echo 'source ~/anaconda3/bin/activate python2' >> ~/.bash_profile

# fix homebrew:
sudo chown -R $(whoami) /usr/local/include
sudo chown -R $(whoami) /usr/local/lib
sudo chown -R $(whoami) /usr/local/lib/pkgconfig
sudo chown -R $(whoami) /usr/local/share/man/man5
sudo chown -R $(whoami) /usr/local/share/man/mann


brew tap homebrew/science
brew install samtools
brew install bcftools

# install Mail ActOn
# install [checkVCF.py](https://github.com/zhanxw/checkVCF)
#    - symlink the script to ~/bin
#    - requires ~/Documents/src/

pip install grip
# this doesn't work no matter how I try to install it: pip install bx-python

