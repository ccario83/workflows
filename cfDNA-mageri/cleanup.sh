#!/bin/bash
echo -e "This will delete all nextflow work files \033[1;37mAND POTENTIALLY OTHER *.log, *.txt, and *.html FILES\033[0m in this and the work/ directories."
echo "This should be OK if you've not modified the workflow directory."
read -r -p "Are you sure? [y/N] " response
case $response in
    [yY][eE][sS]|[yY]) 
        rm -rf work
        rm -f  .nextflow.*
	rm -rf .nextflow 
        rm -f  *timeline.html*
        rm -f  *trace.txt*
        ;;
    *)
        ;;
esac





