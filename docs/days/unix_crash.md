The goal of this crash-course is to get you acquainted with some the basic concepts and commands of the UNIX command line.

Mastering it is not something done in a single day, and so the goal is not to make you into sysadmins, but to give you enough groundings so that,
armed with a good cheatsheet and with a bit of patience, you should be able to run the commands needed to perform a differential expression or enrichment analysis and then be able to interact with your result files.

This crash course is made to be followed along either on a distant Rstudio server (the teacher will tell you how to connect there)


**You will learn to :**

 * navigate the filesystem using common unix commands
 * look at text files and manipulate them 
 * use command line tools options
 * redirect command line tools output to files


## start a terminal in Rstudio

![a view of Rstudio's terminal](../assets/images/rstudio_terminal.png)

## filesystem navigation

type commands and then press enter to execute them

### pwd

`pwd` = **p**rint **w**orking **d**irectory : tells you were the ternminal is in the filesystem  
```sh
pwd
```

### ls

`ls` = **l**i**s**t : list the content of the working directory
```sh
ls
```

**Giving arguments to command - 1 :** some commands take arguments (very much like function in R).

Arguments are given separated by **spaces**.


For example, `ls` can take an argument corresponding the name of a directory, in which case it will list the content of that directory rather than the current working directory:

```sh
ls R_crash_course_data
```



### cd

`cd` = **c**hange **d**irectory : moves the current working directory to the folder given as argument
```sh
cd R_crash_course_data

pwd # check the new working directory
```

.. is a shortcut for "parent directory"
```sh
cd .. 
```

. is a shortcut for "current directory"
```sh
cd . # we move to the current directory, so we stay in place
```

cd without argument will take you back to your home directory
```sh
cd
```

You can move through multiple folders in one go
```sh
cd R_crash_course_data/data/
```


**Giving arguments to command - 2 :**

Some arguments are given by name, in which case they start with `-` (one letter argument) or `--` (more than one letter argument).

Also some arguments require a value, but some are just "flags" to merely need to be added to the command line to take effect:

```sh
ls --color=auto

ls -l
```

One can combine multiple one-letter argument with the same `-`

```R
ls -l -h 

##equivalent to
ls -lh 
```

So in the end you may have a command line looking like 

```sh
ls -lh --color=auto 
```

Its output looks like:
```
total 5.3M
-rw-rw-r-- 1 wandrille wandrille 2.7M Aug 15 13:22 diamonds.csv
-rw-rw-r-- 1 wandrille wandrille 2.7M Aug 15 13:22 diamonds.tsv
-rw-rw-r-- 1 wandrille wandrille 1.4K Aug 15 13:22 virginica.csv
```

* The first line tells you to total size of the files in the folder (but it does not account for sub-folders)

| file permissions | # of links | owner     | group     | size | last modification  | filename     |
|------------------|-----------------|-----------|-----------|------|--------------|--------------|
| -rw-rw-r--       | 1               | wandrille | wandrille | 2.7M | Aug   15 13:22 | diamonds.csv |


File permission are **r**ead **w**rite e**x**ecute

| directory? | owner permissions | group permissions | others permission|
|---|-----|-----|-----|
| - | rw- | rw- | r-- |
| not a directory | can read and write|can read and write|can read only|





## moving files around


First let's move back to our home folder
```sh
cd
```

### mkdir

`mkdir` = **m**a**k**e **dir**ectory : creates a directory
```sh
mkdir unix_exercise

ls -lh # to check that the directory was created
```

### cp

`cp` = **c**o**p**y : copies a file to a given location
```sh
cp R_crash_course_data/iris.csv .  # copying the iris.csv file to the current directory

ls -lh --color=auto  # checking the result
```


Entire folder can be recursively copied, but you need to add the `-r` flag
```sh
 cp R_crash_course_data/data/ unix_exercise/ #trying to copy the data/ folder to unix_exercise/
```
gives the error:
```
cp: -r not specified; omitting directory 'R_crash_course_data/data/'
```

```sh
cp -r R_crash_course_data/data/ unix_exercise/

ls -Rlh --color=auto  unix_exercise ## checking the result; -R makes a recursive ls
```

### mv

`mv` = **m**o**v**e : moves file to a given location
```sh
ls # check that iris.csv is here

mv iris.csv unix_exercise # move iris.csv to the unix_exercise/ folder

ls # check that iris.csv is not here anymore
ls unix_exercise # check that iris.csv is in unix_exercise/
```

it can also be used to change the name of a file or folder
```sh
cd unix_exercise

mv iris.csv IRIS.csv
ls
```

### rm

`rm` = **r**e**m**ove : removes files or folder 

!!! warning
    there is no recycle bin, or trash bin, or whatever.
    **Any deletion is final.**


!!! warning
    I repeat: **any deletion is final.**

let's create a copy of IRIS.csv to delete it:
```sh
cp IRIS.csv iris_copy.csv
ls 
```

```sh
rm IRIS.csv
```

folder's can be deleted with the -r option
```
rm -r data
```

### wildcard character : the * of the show

`*` is the wildcard character, which can subtitute itself to any number of characters.

we use it to craft powerful **regular expressions** that let us apply commands to selected files based on their names:

```sh
# copy all files in ../R_crash_course_data/data/ to the current working directory
cp ../R_crash_course_data/data/* .

ls
```
`diamonds.csv  diamonds.tsv  iris_copy.csv  virginica.csv`

here the `*` substituted itself for all possibles file names.

Where it shine is that it can be combined with other characters:

```
ls -l *.csv ## all csv files

ls -l diamonds.* ## all files starting with diamonds
```


## exercise : basic file system manipulation

 1. change your working directory back to your home directory
 2. create a folder named `unix_exercise2`
 3. look at the files in a folder named `/data/DATA/mouseMT`. What is the size of each file? Do you have permission to write there?
 4. copy all sample_a files from `/data/DATA/mouseMT` to `unix_exercise2`
 5. look at the files you just copied: do you have write permission with them now?



## looking at file content

Let's move to the `unix_exercise` folder
```sh
cd # to your home directory first
cd unix_exercise 
```

### head / tail

`head` shows you the first 10 lines of a file
```sh
head diamonds.csv
```
```
"carat","cut","color","clarity","depth","table","price","x","y","z"
0.23,"Ideal","E","SI2",61.5,55,326,3.95,3.98,2.43
0.21,"Premium","E","SI1",59.8,61,326,3.89,3.84,2.31
0.23,"Good","E","VS1",56.9,65,327,4.05,4.07,2.31
0.29,"Premium","I","VS2",62.4,58,334,4.2,4.23,2.63
0.31,"Good","J","SI2",63.3,58,335,4.34,4.35,2.75
0.24,"Very Good","J","VVS2",62.8,57,336,3.94,3.96,2.48
0.24,"Very Good","I","VVS1",62.3,57,336,3.95,3.98,2.47
0.26,"Very Good","H","SI1",61.9,55,337,4.07,4.11,2.53
0.22,"Fair","E","VS2",65.1,61,337,3.87,3.78,2.49
```

with `-n` you can specify the number of lines to show
```sh
head -n 3 diamonds.csv
```
```
"carat","cut","color","clarity","depth","table","price","x","y","z"
0.23,"Ideal","E","SI2",61.5,55,326,3.95,3.98,2.43
0.21,"Premium","E","SI1",59.8,61,326,3.89,3.84,2.31
```

`tail` shows you the last lines of a file:
```sh
tail -n 3 diamonds.csv
```
```
0.7,"Very Good","D","SI1",62.8,60,2757,5.66,5.68,3.56
0.86,"Premium","H","SI2",61,58,2757,6.15,6.12,3.74
0.75,"Ideal","D","SI2",62.2,55,2757,5.83,5.87,3.64
```

### more 

`more` prints a file content in the terminal.
When the file is large it prints a page, and you can press
 * `enter` : print one more line
 * `space` : print one more page
 * `q` : stop printing


```sh
more diamonds.csv
```


### wc

wc = **w**ord **c**ount : counts the number of lines, words, and characters in a file

```
wc diamonds.csv
```
result: `
  53941   66023 2772143 diamonds.csv
`

We use it most often with option `-l` to get just the number of lines:
```
wc -l iris_copy.csv
```
result: 
`151 iris_copy.csv`


## exercise : lookoing at files

 1. look at the first and last lines of `/data/DATA/mouseMT/sample_a1.fastq`
 2. we have the logs of a Quality Control tool in file `/data/Solutions/Liu2015/010_l_fastqc_Liu2015.o`. Look at its content with more or less : does it look like all files where processed OK?
 3. how many lines are in each of the fastq files `/data/DATA/mouseMT/` (hint: use * to make your life easier)


## get --help

Of course, you cannot guess all commands arguments. 

You can consult the **man**ual page of some command with `man`, but not all tools have one...

Instead, in general you can get a list of arguments with `--help` (or sometimes `-h`)

```
ls --help
```



## advanced but useful stuff

### output redirection : >

`>` let's you make it so that all standard output (but not error messages) of a command go to a file you specify

This is very useful when you call a tool whose logs you would like to keep (the logs of many bioinformatics tools can be input files for other tools).


```
wc -l /data/DATA/mouseMT/*.fastq > mouseMT.file_sizes.txt

more mouseMT.file_sizes.txt
```

### emergency exit : Ctrl+C

Sometimes you will launch a job, which may take a while and you will realize, with horror, that you made a mistake and this is not what you wanted to do!

`Ctrl+C` stops the execution of whatever command is curently running.

For example, the following command counts number of lines in big files, which takes a while.
Let it run for a couple of seconds and then stop it with `Ctrl+C`
```sh
wc -l /data/DATA/Liu2015/*
```

### background job : &

Ooften you want to launch a job which you know will take a while, but you'd like to be able to continue working while it runs. 

By adding `&` at the end of the command line the job will be executed in the background. 


```sh
wc -l /data/DATA/Liu2015/* > Liu2015.file_sizes.txt &
```

You can use `top` or `ps -u` to see your job running (with `top`, type `q` to exit).

When the command is finished you will get in your terminal something like:

```sh
[1]+  Done                    wc -l /data/DATA/Liu2015/* > Liu2015.file_sizes.txt
```

### variables

As in R, is it possible to store data in a variable.

It is mostly useful to help us manipulate file names more easily.

They are declared with `=` , but you have to be careful: use no spaces 

```
INPUT_FOLDER=/data/DATA/mouseMT
```

then you can use them with ${variable_name}:


```sh
ls ${INPUT_FOLDER}/*

ls ${INPUT_FOLDER}/sample_a*
```

!!! note 
    You will notice most of my bash variable are UPPERCASE. It is not mandatory but I find it a useful convention.

!!! note 
    instead of `${INPUT_FOLDER}`, we could write $INPUT_FOLDER but it is a bit more prone to mixup when integrating the variable content within other character, so beware.



### loop

Loops let us "easily" deploy the same command while changing 1 variable.

This is very useful to deploy a bioinformatic tool on many input files in a tidy way.

For example, imagine we want to deploy a command separately on the mouseMT fatsq files of sample a, and of sample b, it could look something like:

```sh
# I use wc -l as my command here

wc -l /data/DATA/mouseMT/sample_a*.fastq > mouseMT.a.file_sizes.txt

wc -l /data/DATA/mouseMT/sample_b*.fastq > mouseMT.b.file_sizes.txt
```

That works, but imagine that you have more than a handful of options and it becomes tedious.

Additionnally, speaking from experience it is easy to make small typos that will riun your day here (eg, adapting the input file names, but not the output file name )

A for loop in bash looks like:

```sh
for VAR in a b
do # done starts the loop block

 ## some command, which will be repeated, ${VAR} will have a different value each time
 wc -l /data/DATA/mouseMT/sample_${VAR}*.fastq > mouseMT.${VAR}.file_sizes.txt

done # done finishes the loop block
```

