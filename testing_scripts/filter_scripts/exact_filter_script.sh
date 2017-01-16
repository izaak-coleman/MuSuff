idx="_hg19SMuFinIdx"
filter=".fil"
fastq=".fastq"
tar_gz=".tar.gz"


for file in `ls $1`
do
# decompress data file
  printf "Decompressing $file \n"
  tar -zxvf $file

# run aligner 
  base_name=`echo ${file:0:(-13)}`
  gunzipped_file=`echo ${file:0:(-7)}`

  printf "Exact match filtering $gunzipped_file\n"
  printf "/homes/ic711/MuSuff/filter/exact_filter/exact_filter /data/ic711/insilico_data/genome_indexes/hg19.fa $gunzipped_file\n"
  /homes/ic711/MuSuff/filter/exact_filter/exact_filter /data/ic711/insilico_data/genome_indexes/hg19.fa $gunzipped_file
  printf "moving  $gunzipped_file$filter to $base_name$idx$fastq$filter \n"
  mv $gunzipped_file$filter $base_name$idx$fastq$filter

# recompress original file
  printf "Recompressing $gunzipped_file \n"
  tar -zcvf  $gunzipped_file$tar_gz $gunzipped_file # recompress unfiltered file

  printf "Recompressing $base_name$idx$fastq$filter \n"
  tar -zcvf $base_name$idx$fastq$filter$tar_gz $base_name$idx$fastq$filter

  printf "\n\n\n"
done

#/homes/ic711/MuSuff/filter/exact_filter/exact_filter /data/ic711/insilico_data/genome_indexes/hg19rec.fa $gunzipped_file
