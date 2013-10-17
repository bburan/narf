#!/bin/bash

# Under linux, you'll need to make sure all of these things are installed
sudo apt-get install subversion git xterm

mkdir matlab
cd matlab

# To get a copy of baphy
svn checkout --username slee --password feartheferret svn+ssh://slee@bhangra.isr.umd.edu/home/svn/baphy

# To get a copy of the 
scp capybara:~/matlab/baphy/Config/BaphyConfigPath.m ~/matlab/baphy/Config/
# Or just edit the file manually so it points to 'lbhb' and uses 'hyrax.ohsu.edu'

# To get a copy of narf, assuming that you have Otter's share 
git clone /auto/data/code/narf-bare/ narf
cp /auto/data/code/

# Then make sure your SSH keys are on the destination machine.

# Finally, make sure that your .bash_login is copied over too