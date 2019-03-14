#! /bin/bash

actual=/tigress/tcomi/aclark4_temp/results/analysis_test/
expected=/tigress/tcomi/aclark4_temp/results/analysisp4e2/
echo starting comarison of $(basename $actual) to $(basename $expected)

for file in $(ls ${expected}*_quality.txt); do
    act=$(echo $file | sed 's/p4e2/_test/g')
    cmp $act $file && echo $file passed! #|| exit
done

echo starting on .fa.gz...
for file in $(ls ${expected}regions/*.fa.gz); do
    act=$(echo $file | sed 's/p4e2/_test/g')
    cmp <(zcat $file) <$act || echo $file failed
done
echo done!
