#/bin/sh
## collect final results into a single directory
find . -name 'SampleID' | \
    while read f; do 
        echo $(cat $f) $(dirname $(dirname $f)) $(ls -lh $(dirname $f)/../*.tar)
    done | \
        sort > result_list.txt
cat result_list.txt | \
    while read SampleID dir rest; do 
        echo $SampleID $dir
        mkdir final/$SampleID
        mv $dir/out $dir/*.tar $dir/*.html final/$SampleID
    done
