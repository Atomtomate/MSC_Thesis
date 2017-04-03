for i in *_GLoc_MF*; do
    awk '$1 > 0' $i > $i.2
done
