#!/bin/bash

parselog_csv() {
# Recursively find all *.log files in the current working directory
# for each, save file path,final step reached, and ns/day to ./status.txt for quick human readable checking,
# and save the same data in a compressed csv format to ./status.csv for compiling datasets
# good for checking the status of everything in a project really quickly, or for compiling benchmarking data
# If benchmarking, I suggest saving relevant metadata in your *.log filenames

  now=$(date '+%Y-%m-%d %H:%M:%S')

  status_file="status.txt"
  csv="status.csv"

  if [ -f $csv ]; then
    echo "Warning: output file $csv already exists"
    echo "appending data to the end of $csv"
    echo "You may wish to check your data"
    echo "DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY" >> $csv
  else
    echo "File,Date,Time,Step,Time(ps),(ns/day),(hour/ns)" >> $csv
  fi
  touch $status_file

  echo "================================" >> $status_file
  date >> $status_file
  echo "================================" >> $status_file


#  for file in $(find . -type f -name '*wet_EQ.log' | sort -k 1) ;do

  for file in $(find . -type f -name '*.log' -not -name '[0-9][0-9]*.log' | sort -k 1) ;do # excludes log files with names like "03*.log" which are usually energy minimisation or eq

    echo "Processing $file"

    # Add the file name to the status file
    echo $file >> $status_file
    echo >> $status_file

    # Get the last line that matches the pattern "Step           Time"
    step_time_line=$(grep -nP 'Step\s+Time' $file | tail -1 | cut -d ":" -f1)
    #echo $step_time_line
    # Add the found line, and the next line, to the status file
    #echo $(tail -2 $file | head -1) >> $status_file
    if [ -n "$step_time_line" ]; then
      foo=$(eval sed -n "$((${step_time_line}+1))p" $file)
      sed -n "${step_time_line}p" $file >> $status_file
      sed -n "$((${step_time_line}+1))p" $file  >> $status_file
      echo >> $status_file
    else
      foo=" NaN   NaN"  
    fi

    # Get the last line that matches the pattern "Performance:"
    #echo $perf_line
    perf_line=$(grep -nP 'Performance:' $file | tail -1| cut -d ":" -f1)
    # Add the found line, and the preceeding line, to the status file
    if [ -n "$perf_line" ]; then
      bar=$(echo $(eval sed -n "${perf_line}p" $file ) | sed 's/Performance://')
      sed -n "$((${perf_line}-1))p" $file >> $status_file
      sed -n "${perf_line}p" $file >> $status_file
    else
      bar=" NaN NaN"
    fi
    echo >> $status_file
    echo >> $status_file
    csv_line=$( echo "$file $now $foo $bar" | sed 's/ \{1,\}/,/g')
    echo $csv_line >> $csv
  done
}
