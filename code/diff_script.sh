while read SAMPLE; do

	cd ${SAMPLE}_alignments/
	echo "Decompressing ${SAMPLE}:"
	gzip -d *.bed.gz
	echo "Running diffReps for ${SAMPLE}:"
	diffReps.pl -tr OD*.bed -co YD*.bed -chrlen hg19_chr.txt -re diff_${SAMPLE}.nb.txt -me nb --nproc 5
	echo "Done"
	echo "Compressing back:"
	gzip *.bed
	cd ../
done < samples.txt


