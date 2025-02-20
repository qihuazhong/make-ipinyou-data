
forked from [wnzhang/make-ipinyou-data](https://github.com/wnzhang/make-ipinyou-data)

Changes:
- Adapted the code to python3
- Only process season 3 data
- Minor correction the data description

## Troubleshooting:

External solvers are needed [(link)](https://yalmip.github.io/allsolvers/) to solve the SDP on YALMIP. Here are steps to install [MOSEK](https://www.mosek.com) (this solver is academic-free but licences are required).
### Step 1: download the solver
Please use the more preferred way for General setup:

link: https://www.mosek.com/downloads/
### Step 2: license
Apply for [Personal Academic License](https://www.mosek.com/products/academic-licenses/). An email will be sent to your mailbox with license attached within 1 minute.

The license should be put inside a folder called "mosek" under the user's home directory. Specifically,

Windows:

```c:\users\_userid_\mosek\mosek.lic```

Unix/OS X:

```/home/_userid_/mosek/mosek.lic```

Where _userid_ is your User ID on the computer.
### Step 3: install the solver 
unpack it into a chosen directory, this path recommended: 

```/home/_userid_/mosek/~```

Run the command

```python <MSKHOME>/mosek/9.1/tools/platform/osx64x86/bin/install.py```
### Step 4: toolbox installation on Matlab 
add path:

```addpath ~/mosek/9.1/toolbox/r2015a```

set the solver:

```ops = sdpsettings('solver','mosek','verbose',1,'debug',1)```

make-ipinyou-data
=================

This project is to formalise the iPinYou RTB data into a standard format for further researches.

### Step 0

The raw data of iPinYou (`ipinyou.contest.dataset.zip`) can be downloaded from [Dropbox](https://www.dropbox.com/s/txz0ms0axqf7jrl/ipinyou.contest.dataset.7z?dl=0).

Unzip it and get the folder `ipinyou.contest.dataset`.

### Step 1
Update the soft link for the folder `ipinyou.contest.dataset` in `original-data`. 
```
weinan@ZHANG:~/Project/make-ipinyou-data/original-data$ ln -sfn ~/Data/ipinyou.contest.dataset ipinyou.contest.dataset
```
Under `make-ipinyou-data/original-data/ipinyou.contest.dataset` there should be the original dataset files like this:
```
weinan@ZHANG:~/Project/make-ipinyou-data/original-data/ipinyou.contest.dataset$ ls
algo.submission.demo.tar.bz2  README         testing2nd   training3rd
city.cn.txt                   region.cn.txt  testing3rd   user.profile.tags.cn.txt
city.en.txt                   region.en.txt  training1st  user.profile.tags.en.txt
files.md5                     testing1st     training2nd
```
You do not need to further unzip the packages in the subfolders.

### Step 2
Under `make-ipinyou-data` folder, just run `make all`.

> :warning: This version of makefile only process season 3 data (2261, 2259, 2821 and 2993)

After the program finished, the total size of the folder will be 14G. The files under `make-ipinyou-data` should be like this:
```
weinan@ZHANG:~/Project/make-ipinyou-data$ ls
1458  2261  2997  3386  3476  LICENSE   mkyzxdata.sh   python     schema.txt
2259  2821  3358  3427  all   Makefile  original-data  README.md
```
Normally, we only do experiment for each campaign (e.g. `1458`). `all` is just the merge of all the campaigns. You can delete `all` if you think it is unuseful in your experiment.

### Use of the data
We use campaign 1458 as example here.
```
weinan@ZHANG:~/Project/make-ipinyou-data/1458$ ls
featindex.txt  test.log.txt  test.yzx.txt  train.log.txt  train.yzx.txt
```
* `train.log.txt` and `test.log.txt` are the formalised string data for each row (record) in train and test. The first column is whether the user click the ad or not. The 14th column is the winning price for this auction.
* `featindex.txt`maps the features to their indexes. For example, `8:115.45.195.*	29` means that the 8th column in `train.log.txt` with the string `115.45.195.*` maps to feature index `29`.
* `train.yzx.txt` and `test.yzx.txt` are the mapped vector data for `train.log.txt` and `test.log.txt`. The format is y:click, z:wining_price, and x:features. Such data is in the standard form as introduced in [iPinYou Benchmarking](http://arxiv.org/abs/1407.7073).


For any questions, please report the issues or contact [Weinan Zhang](http://www0.cs.ucl.ac.uk/staff/w.zhang/).
