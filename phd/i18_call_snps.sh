#!/bin/bash

while getopts 'd:' OPTION; do
  case "$OPTION" in
    d)
      d="$OPTARG"
      echo "The project directory is $OPTARG"
      ;;
    ?)
      echo "script usage: $(basename $0) [-t string]" >&2
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"

### Make project & sub directories as required:
echo "Making project directory if it does not already exist..."
if [ ! -d "$d" ]; then
  mkdir $d
  echo "${d} Project directory created."
else
  echo "${d} Project directory already exists."
fi

echo "Changing directories into project directory..."
cd $d

echo "Making subdirectories if they do not already exist..."
if [ ! -d "${d}/raw_snippy_output/" ]; then
  mkdir "${d}/raw_snippy_output/"
  echo "${d}/raw_snippy_output/ dir created"
else
  echo "raw_snippy_output/ dir exists"
fi

if [ ! -d "${d}/snippy_core/" ]; then
  mkdir "${d}/snippy_core/"
  echo "${d}/snippy_core/ dir created"
else
  echo "snippy_core/ dir exists"
fi

if [ ! -d "${d}/iqtree/" ]; then
  mkdir "${d}/iqtree/"
  echo "${d}/iqtree/ dir created"
else
  echo "iqtree/ dir exists"
fi

if [ ! -d "${d}/cfml/" ]; then
  mkdir "${d}/cfml/"
  echo "${d}/cfml/ dir created"
else
  echo "cfml/ dir exists"
fi

### Continue script...
