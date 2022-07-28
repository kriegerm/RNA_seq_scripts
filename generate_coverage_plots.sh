#Generate Coverage Curves
#change the locaiton of "coverage_curves_PE.pl" as needed.

for file in ./*.bam
do
echo $file
filename=$(basename "$file") 
filename="${filename%.*}"
echo "converting to .sam format..."
samtools view -h $file > $filename.sam
echo "generating coverage plot..."
perl ./coverage_curves_PE.pl $filename.sam > $filename.txt
rm $filename.sam
done
