
# read working directory

pwd > pwd.tmp
read WD < pwd.tmp
rm pwd.tmp

# go in the dir with main.tex and compile

cd ../ 

latex main.tex > tmp

# go back

cd $WD

