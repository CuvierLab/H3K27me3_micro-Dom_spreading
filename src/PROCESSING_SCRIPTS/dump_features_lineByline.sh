# This script aims to do the APA computation to get the APA matrix of a given set of proteins.
# For each protein, the APA is the agregation of a contact score in the form of a matrix of XBP x XBP at a given resolution
# center on the protein of interest.
# We will use then those matrices to compare the signal on center/controls(border which are not the center of the matrix)
# with a Fisher enrichment test

# We have to give a 2d annotation file that represent the size of matrix we want to observe contacts in :
#chr1   x1         x2         chr2   y1         y2         color     comment
#chrX   85000000   89000000   chrX   85000000   89000000   0,255,0   My green region
#chrX   90000000   99100000   chrX   90000000   99100000   0,0,255   My blue region


# juicebox dump <observed/oe> <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize> [outfile]
#          dump <norm/expected> <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr> <BP/FRAG> <binsize> [outfile]
#          dump <loops/domains> <hicFile URL> [outfile]

# GET OPTION
while getopts "f:h:o:r:" option
do
case $option in
    f)
        file=$OPTARG
        ;;
    h)
        hic=$OPTARG
        ;;
    r)
        res=$OPTARG
        ;;
    o)
        outdir=$OPTARG
        ;;
    \?)
        echo "-f emplacement of 2D annotations eg:/my/file.txt"
        echo "-h .hic of reference"
        echo "-r resolution of exctracted matrice (/!\ not all resolution available for a given .hic)"
        echo "-o where the output file should be written eg:/my/idx_ref/directory"
        exit 1
        ;;
    h)
				echo "-f emplacement of 2D annotations eg:/my/file.txt"
				echo "-h .hic of reference"
        echo "-r resolution of exctracted matrice (/!\ not all resolution available for a given .hic)"
				echo "-o where the output file should be written eg:/my/idx_ref/directory"
        ;;
esac
done

#Make output directory
# Select juicer path
juicer=/media/alexandre/Data/Softwares/juicer_tools_0.7.5.jar

# Run the for loop at different resolutions

b_file=`basename $file`
out=${b_file%*.bed}
mkdir -p  ${outdir}

echo ${out}
declare -i cnt=1
while read chr1 s1 s2 chr2 e1 e2
do
  echo "$chr1 $chr2"
	java -jar ${juicer} dump oe KR $hic ${chr1}:${s1}:${s2} ${chr2}:${e1}:${e2} BP ${res} ${outdir}/${out}_${cnt}_${chr1}.txt
  ((cnt++))
  echo ${outdir}/${out}.txt
done < $file
