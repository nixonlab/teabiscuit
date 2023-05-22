# Module 0: Setup

### Login information

Go to Rstudio instance:

[http://nixonlab-teabiscuit.hopto.org/](http://nixonlab-teabiscuit.hopto.org/)

Enter username and password.


### Install conda

```
cd $HOME
ln -s /fsx/users/$(whoami) myfsx
bash /efs/software/Miniconda3-py39_4.10.3-Linux-x86_64.sh -b -p $HOME/miniconda3
$HOME/miniconda3/bin/conda init
```

Now you have to log out and log back in.

### Install mamba

```
conda install -n base -y -c conda-forge mamba
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Create teabiscuit environment:

```
git clone https://github.com/nixonlab/teabiscuit.git
mamba env create -f teabiscuit/00-setup/teabiscuit.yml
conda activate teabiscuit
pip install git+https://github.com/mlbendall/telescope.git
telescope --version
```


