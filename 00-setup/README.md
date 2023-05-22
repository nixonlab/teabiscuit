# Module 0: Setup

### Login information

Go to Rstudio instance:

[http://nixonlab-teabiscuit.hopto.org/](http://nixonlab-teabiscuit.hopto.org/)

Enter username and password.


### Conda first-time setup

```
cd $HOME
ln -s /fsx/users/$(whoami) myfsx
bash /efs/software/Miniconda3-py39_4.10.3-Linux-x86_64.sh -b -p $HOME/miniconda3
$HOME/miniconda3/bin/conda init
```

Now you have to log out and log back in.

#### Install mamba and set up channels

```
conda install -n base -y -c conda-forge mamba
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Useful conda commands

[Conda cheatsheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf)

List your environments:

```
conda info --envs
```

Delete an environment:

```
conda remove -n [envname] --all
```

Use mamba for installing and creating

```
mamba install numpy
mamba create -n newenv python numpy pandas
```


### Create teabiscuit environment:

```
git clone https://github.com/nixonlab/teabiscuit.git
conda env create -f teabiscuit/00-setup/teabiscuit.yml
conda activate teabiscuit
telescope --version
```


