#!/bin/bash

# 逐行读取文件，对每一行进行列的重新排列
while IFS= read -r line; do
    echo "$line" | awk '{print $2, $3, $1}' >> tmp
done < output.txt
mv tmp output.txt

array=$(sed -n '4p' input)
IFS=',' read -r -a array <<< "$array"

readarray lines < output.txt

for i in "${!array[@]}"; do
  prefix="${array[$i]}"
  echo -n "${prefix} ${lines[$i]}" >> temp
done
mv temp output.txt

