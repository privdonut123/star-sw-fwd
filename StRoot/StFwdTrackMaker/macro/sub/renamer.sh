#!/usr/bin/env bash

for file in *.picoDst.root; do
    [[ "$file" =~ _J[^/]+\.picoDst\.root$ ]] && continue
    base="${file%.picoDst.root}"
    echo mv "$file" "${base}_J${1}.picoDst.root"
    mv "$file" "${base}_J${1}.picoDst.root"
done

for file in *.FwdFitQA.root; do
    [[ "$file" =~ _J[^/]+\.FwdFitQA\.root$ ]] && continue
    base="${file%.FwdFitQA.root}"
    echo mv "$file" "${base}_J${1}.FwdFitQA.root"
    mv "$file" "${base}_J${1}.FwdFitQA.root"
done